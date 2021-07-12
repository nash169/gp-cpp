#include <fstream>
#include <iostream>

#include <mfem.hpp>

#include <utils_cpp/UtilsCpp.hpp>

using namespace std;
using namespace mfem;

int main(int argc, char* argv[])
{
    // Params
    const char* mesh_file = "rsc/mesh/sphere.msh";
    int ser_ref_levels = 1, order = 1, nev = 30;
    double dbc_eig = 1e3;
    bool visualization = 1;

    // Mesh
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

    // Refinement
    // for (int lev = 0; lev < ser_ref_levels; lev++) {
    //     mesh->UniformRefinement();
    // }

    // Elements
    FiniteElementCollection* fec;
    if (order > 0) {
        fec = new H1_FECollection(order, dim);
    }
    else if (mesh->GetNodes()) {
        fec = mesh->GetNodes()->OwnFEC();
    }
    else {
        fec = new H1_FECollection(order = 1, dim);
    }
    FiniteElementSpace* fespace = new FiniteElementSpace(mesh, fec);

    // Mass and Stiffness matrix
    ConstantCoefficient one(1.0);
    Array<int> ess_bdr;
    if (mesh->bdr_attributes.Size()) {
        ess_bdr.SetSize(mesh->bdr_attributes.Max());
        ess_bdr = 1;
    }

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

    // Eigenvalue problem
    SpectraEigenSolver spectra;

    spectra.SetNumModes(nev)
        .SetKrylov(50)
        .SetMaxIter(5000)
        .SetTol(1e-5)
        .SetOperators(*a, *m)
        .Solve();

    Eigen::VectorXd eigenvalues = spectra.GetEigenvalues(nev);

    utils_cpp::FileManager io_manager;
    Eigen::IOFormat precision(10);
    io_manager.setFile("rsc/sol/sphere_eigs").write(eigenvalues.format(precision));

    // Save
    GridFunction x(fespace);

    {
        ostringstream mesh_name, mode_name;
        mesh_name << "rsc/sol/"
                  << "sphere"
                  << ".mesh";

        ofstream mesh_ofs(mesh_name.str().c_str());
        mesh_ofs.precision(8);
        mesh->Print(mesh_ofs);

        for (int i = 0; i < nev; i++) {
            // conver Eigen Vector to MFEM Vector
            Vector eigenvector = VectorConverter<double>::from(spectra.GetEigenvector(i));

            // convert eigenvector from Vector to GridFunction
            x = eigenvector;

            mode_name << "rsc/sol/"
                      << "sphere"
                      << "_mode_" << i;

            ofstream mode_ofs(mode_name.str().c_str());
            mode_ofs.precision(8);
            x.Save(mode_ofs);
            mode_name.str("");
        }
    }

    // Plot
    if (visualization) {
        char vishost[] = "localhost";
        int visport = 19916;
        socketstream mode_sock(vishost, visport);
        mode_sock.precision(8);

        for (int i = 0; i < nev; i++) {
            cout << "Eigenmode " << i + 1 << '/' << nev
                 << ", Lambda = " << eigenvalues[i] << endl;

            // convert eigenvector from HypreParVector to ParGridFunction
            Vector eigenvector = VectorConverter<double>::from(spectra.GetEigenvector(i));
            x = eigenvector;

            mode_sock << "solution\n"
                      << *mesh << x << flush
                      << "window_title 'Eigenmode " << i + 1 << '/' << nev
                      << ", Lambda = " << eigenvalues[i] << "'" << endl;

            char c;
            cout << "press (q)uit or (c)ontinue --> " << flush;
            cin >> c;

            if (c != 'c') {
                break;
            }
        }
        mode_sock.close();
    }

    // Free memory
    delete m;
    delete a;

    delete fespace;
    if (order > 0) {
        delete fec;
    }
    delete mesh;

    return 0;
}
