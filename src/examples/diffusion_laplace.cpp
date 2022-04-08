#include <iostream>
#include <sstream>
#include <string>

#include <gp_manifold/SlepcSolver.hpp>
#include <kernel_lib/Kernel.hpp>
#include <utils_lib/FileManager.hpp>

using namespace gp_manifold;
using namespace kernel_lib;
using namespace utils_lib;

// magnitude of the eigenvalues related to the curvature in the different direction of the embedding
// sphere -> -2.99573
// ellipse -> -4.6052
#define SIGMA -4.6052 // -4.6052 -2.99573 -2.30259 -0.6931 (0.01 0.05 0.1 0.5)

struct ParamsExp {
    struct kernel : public defaults::kernel {
        PARAM_SCALAR(double, sf, 0);
        PARAM_SCALAR(double, sn, -5);
    };

    struct exp_sq : public defaults::exp_sq {
        PARAM_SCALAR(double, l, SIGMA);
    };
};

int main(int argc, char** argv)
{
    std::string mesh_name = (argc > 1) ? argv[1] : "sphere",
                mesh_ext = "msh";

    // Load "sampled" nodes
    FileManager io_manager;
    Eigen::MatrixXd vertices = io_manager.setFile("rsc/meshes/" + mesh_name + "." + mesh_ext).read<Eigen::MatrixXd>("$Nodes", 2, "$EndNodes");

    vertices.block(0, 0, vertices.rows(), 3) = vertices.block(0, 1, vertices.rows(), 3);
    vertices.conservativeResize(vertices.rows(), 3);

    // io_manager.setFile("rsc/samples.csv").write(vertices);

    // Geometric Laplacian
    double eps = 2 * std::pow(std::exp(SIGMA), 2);
    size_t num_samples = vertices.rows(), nn = 50, skip_diag = 1;
    utils::Graph graph;

    Eigen::SparseMatrix<double, Eigen::RowMajor> L = graph.kNearestWeighted(vertices, kernels::SquaredExp<ParamsExp>(), nn, skip_diag);

    std::vector<Eigen::Triplet<double>> tripletList;
    for (size_t i = 0; i < L.rows(); i++)
        tripletList.push_back(Eigen::Triplet<double>(i, i, 1 / L.row(i).sum()));
    Eigen::SparseMatrix<double, Eigen::RowMajor> D(num_samples, num_samples);
    D.setFromTriplets(tripletList.begin(), tripletList.end());

    L = D * L * D;

    tripletList.clear();
    for (size_t i = 0; i < L.rows(); i++)
        tripletList.push_back(Eigen::Triplet<double>(i, i, 1 / L.row(i).sum()));
    D.setFromTriplets(tripletList.begin(), tripletList.end());

    L = (Eigen::MatrixXd::Identity(num_samples, num_samples) - D * L) / eps * 4;
    // std::cout << Eigen::MatrixXd(L).block(0, 0, 10, 10) << std::endl;
    // L = (D - L) / eps / 4;

    // io_manager.setFile("rsc/laplacian.csv").write(Eigen::MatrixXd(L));

    // Create Slepc solver
    int nev = (argc > 2) ? std::stoi(argv[2]) : 10;
    SlepcSolver solver(argc, argv);

    // solver.operators(std::make_unique<PetscMatrix>(L), std::make_unique<PetscMatrix>(D)) // Set opeartors
    solver.operators(std::make_unique<PetscMatrix>(L)) // Set opeartors
                                                       // .problem(EPS_NHEP) // Set problem type
        .target(0.0) // target eigenvalue
        .spectrum(EPS_TARGET_REAL) // Select smallest eigenvalues
        .modes(nev) // Number of requested eigenvalues
        .spectralTransformation(STSINVERT)
        // .iterations(500) // Set max iter
        // .tolerance(1e-14) // Set tolerance
        .solve(); // Solve the problem

    // Number of eigenvalues converged
    std::cout << "Number of converged eigenvalues: " << solver.converged() << std::endl;

    // Save eigenvalues & eigenvector
    Eigen::VectorXd eigs(nev);
    Eigen::MatrixXd vecs(num_samples, nev);

    for (size_t i = 0; i < nev; i++) {
        eigs(i) = solver.eigenvalue(i);
        vecs.col(i) = solver.eigenvector(i);
    }

    io_manager.setFile("rsc/modes/diffusion_" + mesh_name + "_modes.000000").write("modes", vecs);
    io_manager.setFile("rsc/modes/diffusion_" + mesh_name + "_eigs.000000").write("eigs", eigs);

    return 0;
}