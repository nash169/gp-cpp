#include <Eigen/Core>
#include <fstream>
#include <iostream>
#include <mfem.hpp>
#include <utils_cpp/UtilsCpp.hpp>

using namespace std;
using namespace mfem;

int main(int argc, char* argv[])
{
    // Initialize MPI
    int num_procs, myid;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    // Data
    string mesh_name = "sphere",
           mesh_ext = "msh",
           mesh_string = "rsc/meshes/" + mesh_name + "." + mesh_ext;

    // Parse options
    const char* mesh_file = mesh_string.c_str();
    int ser_ref_levels = 2;
    int par_ref_levels = 1;
    int order = 1;
    int nev = 5;
    int seed = 75;
    const char* slepcrc_file = "";

    // Initialize SLEPc. This internally initializes PETSc as well.
    MFEMInitializeSlepc(NULL, NULL, slepcrc_file, NULL);

    // Read the (serial) mesh from the given mesh file on all processors. We
    // can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
    // and volume meshes with the same code.
    Mesh* mesh = new Mesh(mesh_file, 1, 1);
    int dim = mesh->Dimension();

    // Refine the serial mesh on all processors to increase the resolution. In
    // this example we do 'ref_levels' of uniform refinement (2 by default, or
    // specified on the command line with -rs).

    // for (int lev = 0; lev < ser_ref_levels; lev++)
    //    mesh->UniformRefinement();

    // Define a parallel mesh by a partitioning of the serial mesh. Refine
    // this mesh further in parallel to increase the resolution (1 time by
    // default, or specified on the command line with -rp). Once the parallel
    // mesh is defined, the serial mesh can be deleted.
    ParMesh* pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;

    // for (int lev = 0; lev < par_ref_levels; lev++)
    //    pmesh->UniformRefinement();

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
    OperatorHandle Ah(Operator::PETSC_MATAIJ), Mh(Operator::PETSC_MATAIJ);

    a->ParallelAssemble(Ah);
    Ah.Get(pA);
    Ah.SetOperatorOwner(false);

    m->ParallelAssemble(Mh);
    Mh.Get(pM);
    Mh.SetOperatorOwner(false);

    delete a;
    delete m;

    // Set the matrices which define the generalized eigenproblem A x = lambda M x.
    SlepcEigenSolver* slepc = NULL;
    slepc = new SlepcEigenSolver(MPI_COMM_WORLD);
    slepc->SetNumModes(nev);
    slepc->SetWhichEigenpairs(SlepcEigenSolver::TARGET_REAL);
    slepc->SetTarget(0.0);
    slepc->SetSpectralTransformation(SlepcEigenSolver::SHIFT_INVERT);
    slepc->SetOperators(*pA, *pM);

    // Compute the eigenmodes and extract the array of eigenvalues. Define a
    // parallel grid function to represent each of the eigenmodes returned by
    // the solver.
    Array<double> eigenvalues;
    Eigen::VectorXd eigs(nev);
    slepc->Solve();
    eigenvalues.SetSize(nev);
    for (int i = 0; i < nev; i++) {
        slepc->GetEigenvalue(i, eigenvalues[i]);
        eigs(i) = eigenvalues[i];
    }

    // Save eigenvalues
    utils_cpp::FileManager io_manager;
    io_manager.setFile("rsc/modes/fem_" + mesh_name + "_eigs.000000").write(eigs);

    Vector temp(fespace->GetTrueVSize());
    ParGridFunction x(fespace);

    // Save the refined mesh and the modes in parallel. This output can be
    // viewed later using GLVis: "glvis -np <np> -m mesh -g mode".
    {
        ostringstream mesh_path, mode_name;
        mesh_path << "rsc/modes/fem_" << mesh_name << "_mesh." << setfill('0') << setw(6) << myid;

        ofstream mesh_ofs(mesh_path.str().c_str());
        mesh_ofs.precision(8);
        pmesh->Print(mesh_ofs);

        for (int i = 0; i < nev; i++) {
            slepc->GetEigenvector(i, temp);
            x.Distribute(temp);

            mode_name << "rsc/modes/fem_" << mesh_name << "_mode_" << i << "."
                      << setfill('0') << setw(6) << myid;

            ofstream mode_ofs(mode_name.str().c_str());
            mode_ofs.precision(8);
            x.Save(mode_ofs);
            mode_name.str("");
        }
    }

    // Free the used memory.
    delete slepc;
    delete pA;
    delete pM;
    delete fespace;
    if (order > 0) {
        delete fec;
    }
    delete pmesh;

    // We finalize SLEPc
    MFEMFinalizeSlepc();
    MPI_Finalize();

    return 0;
}