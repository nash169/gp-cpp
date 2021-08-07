#include <iostream>
#include <sstream>
#include <string>

#include <Eigen/Core>

#include <gp_manifold/GaussianProcess.hpp>
#include <gp_manifold/SlepcSolver.hpp>

#include <magnum_dynamics/MagnumApp.hpp>
#include <utils_cpp/UtilsCpp.hpp>

using namespace gp_manifold;
using namespace kernel_lib;

struct ParamsExp {
    struct kernel : public defaults::kernel {
        PARAM_SCALAR(double, sf, 0);
        PARAM_SCALAR(double, sn, -5);
    };

    struct exp_sq : public defaults::exp_sq {
        PARAM_SCALAR(double, l, -2.30259); // -2.99573
    };
};

struct ParamsRiemann {
    struct kernel : public defaults::kernel {
        PARAM_SCALAR(double, sf, 0);
        PARAM_SCALAR(double, sn, -5);
    };

    struct riemann_exp_sq : public defaults::riemann_exp_sq {
        PARAM_SCALAR(double, l, -2.3025);
    };

    struct exp_sq_full : public defaults::exp_sq_full {
    };

    struct gaussian : public defaults::gaussian {
    };
};

int main(int argc, char** argv)
{
    std::string mesh_name = "sphere", mesh_ext = "msh";

    // Load "sampled" nodes
    utils_cpp::FileManager io_manager;
    Eigen::MatrixXd vertices = io_manager.setFile("rsc/meshes/" + mesh_name + "." + mesh_ext).read<Eigen::MatrixXd>("$Nodes", 2, "$EndNodes"),
                    indices = io_manager.read<Eigen::MatrixXd>("$Elements", 2, "$EndElements").array() - 1;

    vertices.block(0, 0, vertices.rows(), 3) = vertices.block(0, 1, vertices.rows(), 3);
    vertices.conservativeResize(vertices.rows(), 3);

    // Load ground truth, target and relative nodes
    Eigen::MatrixXd nodes = io_manager.setFile("rsc/truth/" + mesh_name + "_vertices.csv").read<Eigen::MatrixXd>(),
                    faces = io_manager.setFile("rsc/truth/" + mesh_name + "_faces.csv").read<Eigen::MatrixXd>(),
                    reference = io_manager.setFile("rsc/truth/" + mesh_name + "_reference.csv").read<Eigen::MatrixXd>();

    Eigen::VectorXd ground_truth = io_manager.setFile("rsc/truth/" + mesh_name + "_truth.csv").read<Eigen::MatrixXd>(),
                    target = io_manager.setFile("rsc/truth/" + mesh_name + "_target.csv").read<Eigen::MatrixXd>();

    // Geometric Laplacian
    double eps = 2 * std::pow(std::exp(-2.30259), 2);
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

    L = (Eigen::MatrixXd::Identity(num_samples, num_samples) - D * L) / eps / 4;
    // L = (D - L) / eps / 4;

    // Create Slepc solver
    int nev = 5;
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
    for (size_t i = 0; i < nev; i++) {
        eigs(i) = solver.eigenvalue(i);
        io_manager.setFile("rsc/modes/diffusion_" + mesh_name + "_mode_" + std::to_string(i) + ".000000").write("mode", solver.eigenvector(i));
    }

    io_manager.setFile("rsc/modes/diffusion_" + mesh_name + "_eigs.000000").write(eigs);

    return 0;
}