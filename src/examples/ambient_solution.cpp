#include <iostream>
#include <random>
#include <sstream>

#include <gp_manifold/GaussianProcess.hpp>
#include <utils_lib/FileManager.hpp>

#include <Eigen/Dense>

using namespace gp_manifold;
using namespace kernel_lib;
using namespace utils_lib;

struct ParamsExp {
    struct kernel : public defaults::kernel {
        PARAM_SCALAR(double, sf, 0);
        PARAM_SCALAR(double, sn, -5);
    };

    struct exp_sq : public defaults::exp_sq {
        PARAM_SCALAR(double, l, -0.3566);
    };
};

int main(int argc, char** argv)
{
    std::string mesh_name = (argc > 1) ? argv[1] : "sphere";
    FileManager io_manager;

    // Load ground truth, target and relative nodes
    Eigen::MatrixXd nodes = io_manager.setFile("rsc/truth/" + mesh_name + "_vertices.csv").read<Eigen::MatrixXd>();
    Eigen::VectorXd ground_truth = io_manager.setFile("rsc/truth/" + mesh_name + "_truth.csv").read<Eigen::MatrixXd>();

    // Ambient space Gaussian Process
    using Kernel_t = kernels::SquaredExp<ParamsExp>;
    using GP_t = GaussianProcess<ParamsExp, Kernel_t>;
    GP_t gp;

    constexpr size_t NUM_RUN = 10;
    constexpr int RAND_NUMS_TO_GENERATE[] = {25, 50, 75, 100, 125, 150};

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
            gp.setSamples(reference).setTarget(target).update();

            // Evaluation of the GP on all the mesh points
            gp_sol = gp.multiEval(nodes);

            // record mse
            mse(j, i) = (ground_truth - gp_sol).array().pow(2).sum() / nodes.rows();
        }
    }

    // Save GP solution
    io_manager
        .setFile("rsc/solutions/ambient_" + mesh_name + "_gp.csv")
        .write("mse", mse, "sol", gp_sol);

    return 0;
}