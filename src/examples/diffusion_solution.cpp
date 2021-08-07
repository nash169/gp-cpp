#include <iostream>
#include <sstream>

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
        PARAM_SCALAR(double, l, -0.3566);
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
    double eps = 2 * std::pow(std::exp(-0.3566), 2);
    size_t num_samples = vertices.rows(), nn = 50, skip_diag = 1;
    utils::Graph graph;

    Eigen::SparseMatrix<double, Eigen::RowMajor> L = graph.kNearestWeighted(vertices, kernels::SquaredExp<ParamsExp>(), nn, skip_diag);

    std::vector<Eigen::Triplet<double>> tripletList;
    for (size_t i = 0; i < L.rows(); i++)
        tripletList.push_back(Eigen::Triplet<double>(i, i, 1 / L.row(i).sum()));
    Eigen::SparseMatrix<double, Eigen::RowMajor> D(num_samples, num_samples);
    D.setFromTriplets(tripletList.begin(), tripletList.end());

    L = D * L * D;
    for (size_t i = 0; i < L.rows(); i++)
        tripletList.push_back(Eigen::Triplet<double>(i, i, 1 / L.row(i).sum()));
    D.setFromTriplets(tripletList.begin(), tripletList.end());

    L = (D - L) / eps / 4;

    // Create Slepc solver
    SlepcSolver solver(argc, argv);

    solver.operators(std::make_unique<PetscMatrix>(L), std::make_unique<PetscMatrix>(D)) // Set opeartors
        .problem(EPS_PGNHEP) // Set problem type
        .spectrum(EPS_SMALLEST_MAGNITUDE) // Select smallest eigenvalues
        .modes(5) // Number of requested eigenvalues
        .iterations(100) // Set max iter
        .tolerance(1e-3) // Set tolerance
        .solve(); // Solve the problem

    // // Number of eigenvalues converged
    // std::cout << "Number of converged eigenvalues: " << solver.converged() << std::endl;

    // // Get eigenvector
    // Vec eigenvec = solver.eigenvector(1);

    // PetscPrintf(PETSC_COMM_WORLD, "\nThe eigenvector is:\n\n");
    // VecView(eigenvec, PETSC_VIEWER_STDOUT_WORLD);

    // // Riemann Gaussian Process
    // using Kernel_t = kernels::SquaredExp<ParamsExp>;
    // using Expansion_t = utils::Expansion<ParamsExp, Kernel_t>;
    // using Riemann_t = kernels::RiemannSqExp<ParamsRiemann, Expansion_t>;
    // using RGP_t = GaussianProcess<ParamsRiemann, Riemann_t>;
    // RGP_t rgp;

    // // Set kernel eigen pairs
    // for (size_t i = 1; i < num_eig; i++) {
    //     // Create eigenfunction
    //     Expansion_t f;

    //     // Set manifold sampled points and weights
    //     f.setReference(vertices).setParams(eigenvectors.col(i));

    //     // Add eigen-pair to Riemann kernel
    //     rgp.kernel().addPair(eigenvalues(i), f);
    // }

    // // Set training point and target
    // rgp.setReference(reference).setTarget(target).update();

    // // Evaluation of the GP on all the mesh points
    // Eigen::VectorXd rgp_sol = rgp.multiEval2(nodes);

    // // Save GP solution
    // io_manager.setFile("rsc/solutions/" + mesh_name + "_rgp_diffusion.csv").write(rgp_sol);

    return 0;
}