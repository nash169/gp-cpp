#include <Eigen/Core>
#include <fstream>
#include <iostream>
#include <mfem.hpp>
#include <utils_lib/FileManager.hpp>

using namespace std;
using namespace mfem;
using namespace utils_lib;

int main(int argc, char* argv[])
{
    // Initialize MPI and HYPRE.
    Mpi::Init(argc, argv);
    int num_procs = Mpi::WorldSize();
    int myid = Mpi::WorldRank();
    Hypre::Init();

    // Data
    string mesh_name = (argc > 1) ? argv[1] : "sphere",
           mesh_ext = "msh",
           mesh_string = "rsc/" + mesh_name + "." + mesh_ext;

    // Parse options
    const char* mesh_file = mesh_string.c_str();
    int ser_ref_levels = 2;
    int par_ref_levels = 1;
    int order = 1;
    int nev = (argc > 2) ? std::stoi(argv[2]) : 10;
    int seed = 75;
    bool slu_solver = false;
    bool sp_solver = false;
    bool visualization = 1;
    bool use_slepc = true;
    const char* slepcrc_file = "";
    const char* device_config = "cpu";

    // Initialize SLEPc. This internally initializes PETSc as well.
    MFEMInitializeSlepc(NULL, NULL, slepcrc_file, NULL);

    // Read the (serial) mesh from the given mesh file on all processors. We
    // can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
    // and volume meshes with the same code.
    Mesh* mesh = new Mesh(mesh_file, 1, 1);
    int dim = mesh->Dimension();

    // Define a parallel mesh by a partitioning of the serial mesh. Refine
    // this mesh further in parallel to increase the resolution (1 time by
    // default, or specified on the command line with -rp). Once the parallel
    // mesh is defined, the serial mesh can be deleted.
    ParMesh* pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;

    // Define a parallel finite element space on the parallel mesh. Here we
    // use continuous Lagrange finite elements of the specified order. If
    // order < 1, we instead use an isoparametric/isogeometric space.
    FiniteElementCollection* fec;
    if (order > 0)
        fec = new H1_FECollection(order, dim);
    else if (pmesh->GetNodes())
        fec = pmesh->GetNodes()->OwnFEC();
    else
        fec = new H1_FECollection(order = 1, dim);

    ParFiniteElementSpace* fespace = new ParFiniteElementSpace(pmesh, fec);
    HYPRE_BigInt size = fespace->GlobalTrueVSize();

    if (myid == 0)
        cout << "Number of unknowns: " << size << endl;

    // Set up the parallel bilinear forms a(.,.) and m(.,.) on the finite
    // element space. The first corresponds to the Laplacian operator -Delta,
    // while the second is a simple mass matrix needed on the right hand side
    // of the generalized eigenvalue problem below. The boundary conditions
    // are implemented by elimination with special values on the diagonal to
    // shift the Dirichlet eigenvalues out of the computational range. After
    // serial and parallel assembly we extract the corresponding parallel
    // matrices A and M.
    ConstantCoefficient one(1.0);
    Array<int> ess_bdr;
    if (pmesh->bdr_attributes.Size()) {
        ess_bdr.SetSize(pmesh->bdr_attributes.Max());
        ess_bdr = 1;
    }

    ParBilinearForm* a = new ParBilinearForm(fespace);
    a->AddDomainIntegrator(new DiffusionIntegrator(one));
    if (pmesh->bdr_attributes.Size() == 0) {
        // Add a mass term if the mesh has no boundary, e.g. periodic mesh or
        // closed surface.
        a->AddDomainIntegrator(new MassIntegrator(one));
    }
    a->Assemble();
    a->EliminateEssentialBCDiag(ess_bdr, 1.0);
    a->Finalize();

    ParBilinearForm* m = new ParBilinearForm(fespace);
    m->AddDomainIntegrator(new MassIntegrator(one));
    m->Assemble();
    // shift the eigenvalue corresponding to eliminated dofs to a large value
    m->EliminateEssentialBCDiag(ess_bdr, numeric_limits<double>::min());
    m->Finalize();

    PetscParMatrix *pA = NULL, *pM = NULL;
    HypreParMatrix *A = NULL, *M = NULL;
    Operator::Type tid = !use_slepc ? Operator::Hypre_ParCSR : Operator::PETSC_MATAIJ;
    OperatorHandle Ah(tid), Mh(tid);

    a->ParallelAssemble(Ah);
    if (!use_slepc) {
        Ah.Get(A);
    }
    else {
        Ah.Get(pA);
    }
    Ah.SetOperatorOwner(false);

    m->ParallelAssemble(Mh);
    if (!use_slepc) {
        Mh.Get(M);
    }
    else {
        Mh.Get(pM);
    }
    Mh.SetOperatorOwner(false);

#if defined(MFEM_USE_SUPERLU) || defined(MFEM_USE_STRUMPACK)
    Operator* Arow = NULL;
#ifdef MFEM_USE_SUPERLU
    if (slu_solver) {
        Arow = new SuperLURowLocMatrix(*A);
    }
#endif
#ifdef MFEM_USE_STRUMPACK
    if (sp_solver) {
        Arow = new STRUMPACKRowLocMatrix(*A);
    }
#endif
#endif

    delete a;
    delete m;

    // Set the matrices which define the generalized eigenproblem A x = lambda M x.
    Solver* precond = NULL;
    if (!use_slepc) {
        if (!slu_solver && !sp_solver) {
            HypreBoomerAMG* amg = new HypreBoomerAMG(*A);
            amg->SetPrintLevel(0);
            precond = amg;
        }
        else {
#ifdef MFEM_USE_SUPERLU
            if (slu_solver) {
                SuperLUSolver* superlu = new SuperLUSolver(MPI_COMM_WORLD);
                superlu->SetPrintStatistics(false);
                superlu->SetSymmetricPattern(true);
                superlu->SetColumnPermutation(superlu::PARMETIS);
                superlu->SetOperator(*Arow);
                precond = superlu;
            }
#endif
#ifdef MFEM_USE_STRUMPACK
            if (sp_solver) {
                STRUMPACKSolver* strumpack = new STRUMPACKSolver(argc, argv, MPI_COMM_WORLD);
                strumpack->SetPrintFactorStatistics(true);
                strumpack->SetPrintSolveStatistics(false);
                strumpack->SetKrylovSolver(strumpack::KrylovSolver::DIRECT);
                strumpack->SetReorderingStrategy(strumpack::ReorderingStrategy::METIS);
                strumpack->DisableMatching();
                strumpack->SetOperator(*Arow);
                strumpack->SetFromCommandLine();
                precond = strumpack;
            }
#endif
        }
    }

    HypreLOBPCG* lobpcg = NULL;
    SlepcEigenSolver* slepc = NULL;
    if (!use_slepc) {

        lobpcg = new HypreLOBPCG(MPI_COMM_WORLD);
        lobpcg->SetNumModes(nev);
        lobpcg->SetRandomSeed(seed);
        lobpcg->SetPreconditioner(*precond);
        lobpcg->SetMaxIter(200);
        lobpcg->SetTol(1e-8);
        lobpcg->SetPrecondUsageMode(1);
        lobpcg->SetPrintLevel(1);
        lobpcg->SetMassMatrix(*M);
        lobpcg->SetOperator(*A);
    }
    else {
        slepc = new SlepcEigenSolver(MPI_COMM_WORLD);
        slepc->SetNumModes(nev);
        slepc->SetWhichEigenpairs(SlepcEigenSolver::TARGET_REAL);
        slepc->SetTarget(0.0);
        slepc->SetSpectralTransformation(SlepcEigenSolver::SHIFT_INVERT);
        slepc->SetOperators(*pA, *pM);
    }

    // Compute the eigenmodes and extract the array of eigenvalues. Define a
    // parallel grid function to represent each of the eigenmodes returned by
    // the solver.
    Array<double> eigenvalues;
    if (!use_slepc) {
        lobpcg->Solve();
        lobpcg->GetEigenvalues(eigenvalues);
    }
    else {
        slepc->Solve();
        eigenvalues.SetSize(nev);
        for (int i = 0; i < nev; i++) {
            slepc->GetEigenvalue(i, eigenvalues[i]);
        }
    }

    // Save eigenvalues
    FileManager io_manager;
    Eigen::VectorXd eigs(nev);
    for (int i = 0; i < nev; i++)
        eigs(i) = eigenvalues[i];
    io_manager.setFile("outputs/modes/fem_" + mesh_name + "_eigs.000000").write("eigs", eigs);

    // Save the refined mesh and the modes in parallel. This output can be
    // viewed later using GLVis: "glvis -np <np> -m mesh -g mode".
    Vector temp(fespace->GetTrueVSize());
    ParGridFunction x(fespace);
    Eigen::MatrixXd vecs(size, nev);
    {
        ostringstream mesh_path, mode_name;
        mesh_path << "outputs/modes/fem_" << mesh_name << "_mesh." << setfill('0') << setw(6) << myid;

        ofstream mesh_ofs(mesh_path.str().c_str());
        mesh_ofs.precision(8);
        pmesh->Print(mesh_ofs);

        for (int i = 0; i < nev; i++) {
            if (!use_slepc) {
                x = lobpcg->GetEigenvector(i);
            }
            else {
                slepc->GetEigenvector(i, temp);
                x.Distribute(temp);
            }

            // For now we use directly the vector (works only with one process; check the example to improve this)
            vecs.col(i) = Eigen::Map<Eigen::VectorXd>(x.GetData(), x.Size());
        }
    }
    io_manager.setFile("outputs/modes/fem_" + mesh_name + "_modes.000000").write("modes", vecs);

    // Free the used memory.
    delete lobpcg;
    delete slepc;
    delete precond;
    delete M;
    delete A;
    delete pA;
    delete pM;
#if defined(MFEM_USE_SUPERLU) || defined(MFEM_USE_STRUMPACK)
    delete Arow;
#endif
    delete fespace;
    if (order > 0) {
        delete fec;
    }
    delete pmesh;

    // Finalize SLEPc
    MFEMFinalizeSlepc();

    return 0;
}