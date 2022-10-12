#include <iostream>
#include <random>
#include <sstream>

#include <gp_manifold/GaussianProcess.hpp>
#include <gp_manifold/SlepcSolver.hpp>

#include <utils_lib/FileManager.hpp>

using namespace gp_manifold;
using namespace kernel_lib;
using namespace utils_lib;

struct ParamsExp {
    struct kernel : public defaults::kernel {
        PARAM_SCALAR(double, sf, 0);
        PARAM_SCALAR(double, sn, -5);
    };

    struct exp_sq : public defaults::exp_sq {
        PARAM_SCALAR(double, l, -2.30259);
    };
};

struct ParamsRiemann {
    struct kernel : public defaults::kernel {
        PARAM_SCALAR(double, sf, 0);
        PARAM_SCALAR(double, sn, -5);
    };

    struct riemann_exp_sq : public defaults::riemann_exp_sq {
        PARAM_SCALAR(double, l, -0.6931);
    };
};

int main(int argc, char** argv)
{
    std::string mesh_name = (argc > 1) ? argv[1] : "sphere",
                mesh_ext = "msh";
    FileManager io_manager;

    // Load "sampled" nodes, eigenvectors and eigenvalues for Riemannian kernel construction
    Eigen::MatrixXd vertices = io_manager.setFile("rsc/" + mesh_name + "." + mesh_ext).read<Eigen::MatrixXd>("$Nodes", 2, "$EndNodes"),
                    eigenvectors = io_manager.setFile("outputs/modes/diffusion_" + mesh_name + "_modes.000000").read<Eigen::MatrixXd>("modes", 2);
    Eigen::VectorXd eigenvalues = io_manager.setFile("outputs/modes/diffusion_" + mesh_name + "_eigs.000000").read<Eigen::MatrixXd>("eigs", 2);

    vertices.block(0, 0, vertices.rows(), 3) = vertices.block(0, 1, vertices.rows(), 3);
    vertices.conservativeResize(vertices.rows(), 3);

    // Load ground truth, target and relative nodes
    Eigen::MatrixXd nodes = io_manager.setFile("outputs/truth/" + mesh_name + "_vertices.csv").read<Eigen::MatrixXd>();
    Eigen::VectorXd ground_truth = io_manager.setFile("outputs/truth/" + mesh_name + "_truth.csv").read<Eigen::MatrixXd>();

    // Riemann Gaussian Process
    using Kernel_t = kernels::SquaredExp<ParamsExp>;
    using Expansion_t = utils::Expansion<ParamsExp, Kernel_t>;
    using Riemann_t = kernels::RiemannSqExp<ParamsRiemann, Expansion_t>;
    using RGP_t = GaussianProcess<ParamsRiemann, Riemann_t>;
    RGP_t rgp;

    // Set kernel eigen pairs
    int num_modes = (argc > 2) ? std::stoi(argv[2]) : 10;

    for (size_t i = 1; i < num_modes; i++) {
        // Create eigenfunction
        Expansion_t f;

        // Set manifold sampled points and weights
        f.setSamples(vertices).setWeights(eigenvectors.col(i));

        // Add eigen-pair to Riemann kernel
        rgp.kernel().addPair(eigenvalues(i), f);
    }

    constexpr size_t NUM_RUN = 1;
    constexpr int RAND_NUMS_TO_GENERATE[] = {150}; // {25, 50, 75, 100, 125, 150};

    Eigen::VectorXd gp_sol(nodes.rows());
    Eigen::MatrixXd mse = Eigen::MatrixXd::Zero(NUM_RUN, sizeof(RAND_NUMS_TO_GENERATE) / sizeof(RAND_NUMS_TO_GENERATE[0]));

    for (size_t i = 0; i < sizeof(RAND_NUMS_TO_GENERATE) / sizeof(RAND_NUMS_TO_GENERATE[0]); i++) {
        for (size_t j = 0; j < NUM_RUN; j++) {
            std::cout << RAND_NUMS_TO_GENERATE[i] << " points - trial " << j + 1 << std::endl;

            // Create test set
            Eigen::MatrixXd reference(RAND_NUMS_TO_GENERATE[i], nodes.cols());
            Eigen::VectorXd target(RAND_NUMS_TO_GENERATE[i]);

            std::random_device rd;
            std::default_random_engine eng(rd());
            std::uniform_int_distribution<int> distr(0, nodes.rows() - 1);

            for (size_t k = 0; k < RAND_NUMS_TO_GENERATE[i]; k++) {
                int index = distr(eng);
                reference.row(k) = nodes.row(index);
                target(k) = ground_truth(index);
            }

            // Set training point and target
            rgp.setSamples(reference).setTarget(target).update();

            // Evaluation of the GP on all the mesh points
            gp_sol = rgp.multiEval2(nodes);

            // record mse
            mse(j, i) = (ground_truth - gp_sol).array().pow(2).sum() / nodes.rows();
        }
    }

    // Save GP solution
    io_manager.setFile("outputs/solutions/diffusion_" + mesh_name + "_gp.csv")
        .write("mse", mse, "sol", gp_sol);

    return 0;
}