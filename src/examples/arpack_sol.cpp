#include <fstream>
#include <iostream>
#include <mfem.hpp>

using namespace std;
using namespace mfem;

int main(int argc, char* argv[])
{
    // Mesh default
    std::string name = "torus", extension = ".msh", load_path = "rsc/mesh/", save_path = "rsc/sol/",
                mesh_path = load_path + name + extension;

    const char* mesh_file = mesh_path.c_str();

    // Number of times to refine the mesh (in serial)
    int ref_levels = 3;

    // Shape functions degree
    int order = 1;

    // Number of desired eigenmodes
    int nev = 5;

    // Eigenvalues associated with Dirichlet boundary condition
    double dbc_eig = 1e3;

    // Parse mesh
    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
        "Mesh file to use.");

    // Load the mesh
    Mesh* mesh;
    ifstream imesh(mesh_file);
    if (!imesh) {
        cerr << "\nCan not open mesh file: " << mesh_file << '\n'
             << endl;
        return 2;
    }
    mesh = new Mesh(imesh, 1, 1);
    imesh.close();
    int dim = mesh->Dimension();

    // Refine the mesh
    // for (int lev = 0; lev < ref_levels; lev++)
    //     mesh->UniformRefinement();

    // Shape function order
    FiniteElementCollection* fec;
    fec = new H1_FECollection(order, dim);

    FiniteElementSpace* fespace = new FiniteElementSpace(mesh, fec);

    int size = fespace->GetVSize();

    // Probably some kind of boundary condition for non-compact manifolds
    Array<int> ess_bdr;
    if (mesh->bdr_attributes.Size()) {
        ess_bdr.SetSize(mesh->bdr_attributes.Max());
        ess_bdr = 1;
    }

    // Build stiffness and mass matrix
    ConstantCoefficient one(1.0);

    BilinearForm* a = new BilinearForm(fespace);
    a->AddDomainIntegrator(new DiffusionIntegrator(one));
    if (mesh->bdr_attributes.Size() == 0) {
        // Add a mass term if the mesh has no boundary, e.g. periodic mesh or
        // closed surface.
        a->AddDomainIntegrator(new MassIntegrator(one));
    }
    a->Assemble();
    if (mesh->bdr_attributes.Size() != 0) {
        a->EliminateEssentialBCDiag(ess_bdr, dbc_eig);
    }
    a->Finalize();

    BilinearForm* m = new BilinearForm(fespace);
    m->AddDomainIntegrator(new MassIntegrator(one));
    m->Assemble();
    if (mesh->bdr_attributes.Size() != 0) {
        // shift the eigenvalue corresponding to eliminated dofs to a large value
        m->EliminateEssentialBCDiag(ess_bdr, 1.0);
    }
    m->Finalize();

    // Define and configure the ARPACK eigensolver
    ArPackSym* arpack = new ArPackSym();
    Solver* solver = NULL;

    // Define a simple symmetric Gauss-Seidel preconditioner and use it to solve the system A X = B with PCG.
    GSSmoother M(m->SpMat());
    CGSolver* cg_solver = new CGSolver;
    cg_solver->SetPreconditioner(M);
    cg_solver->SetRelTol(1.0e-12);
    solver = cg_solver;
    solver->SetOperator(m->SpMat());

    // Set arpack
    arpack->SetNumModes(nev);
    arpack->SetMaxIter(400);
    arpack->SetTol(1e-8);
    arpack->SetMode(2);
    arpack->SetPrintLevel(2);

    arpack->SetOperator(*a);
    arpack->SetMassMatrix(*m);
    arpack->SetSolver(*solver);

    // Compute the eigenmodes and extract the array of eigenvalues
    Array<double> eigenvalues;
    arpack->Solve();
    arpack->GetEigenvalues(eigenvalues);

    cout << endl;
    std::ios::fmtflags old_fmt = cout.flags();
    cout.setf(std::ios::scientific);
    std::streamsize old_prec = cout.precision(14);
    for (int i = 0; i < nev; i++) {
        cout << "Eigenvalue lambda   " << eigenvalues[i] << endl;
    }
    cout.precision(old_prec);
    cout.flags(old_fmt);
    cout << endl;

    // Define a grid function to represent each of the eigenmodes returned by the solver
    GridFunction x(fespace);

    // Save the refined mesh and the modes in parallel. This output can be viewed later using GLVis: "glvis -np <np> -m mesh -g mode"
    {
        ostringstream mesh_name, mode_name;
        mesh_name << save_path << name << "_mesh.mesh";

        ofstream mesh_ofs(mesh_name.str().c_str());
        mesh_ofs.precision(8);
        mesh->Print(mesh_ofs);

        for (int i = 0; i < nev; i++) {
            // convert eigenvector from HypreParVector to ParGridFunction
            x = arpack->GetEigenvector(i);

            mode_name << save_path << name << "_mode_" << setfill('0') << setw(2) << i;

            ofstream mode_ofs(mode_name.str().c_str());
            mode_ofs.precision(8);
            x.Save(mode_ofs);
            mode_name.str("");
        }
    }

    // Free the used memory.
    delete arpack;
    delete solver;
    delete m;
    delete a;

    delete fespace;
    if (order > 0) {
        delete fec;
    }
    delete mesh;

    return 0;
}
