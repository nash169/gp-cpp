#include <iostream>
#include <sstream>

#include <Eigen/Core>
#include <gp_manifold/GaussianProcess.hpp>
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

int main(int argc, char** argv)
{
    std::string mesh_name = "sphere";
    utils_cpp::FileManager io_manager;

    // Load ground truth, target and relative nodes
    Eigen::MatrixXd nodes = io_manager.setFile("rsc/truth/" + mesh_name + "_vertices.csv").read<Eigen::MatrixXd>(),
                    faces = io_manager.setFile("rsc/truth/" + mesh_name + "_faces.csv").read<Eigen::MatrixXd>(),
                    reference = io_manager.setFile("rsc/truth/" + mesh_name + "_reference.csv").read<Eigen::MatrixXd>();

    Eigen::VectorXd ground_truth = io_manager.setFile("rsc/truth/" + mesh_name + "_truth.csv").read<Eigen::MatrixXd>(),
                    target = io_manager.setFile("rsc/truth/" + mesh_name + "_target.csv").read<Eigen::MatrixXd>();

    // Ambient space Gaussian Process
    using Kernel_t = kernels::SquaredExp<ParamsExp>;
    using GP_t = GaussianProcess<ParamsExp, Kernel_t>;
    GP_t gp;

    // Set training point and target
    gp.setReference(reference).setTarget(target).update();

    // Evaluation of the GP on all the mesh points
    Eigen::VectorXd gp_sol = gp.multiEval(nodes);

    // Save GP solution
    io_manager.setFile("rsc/solutions/" + mesh_name + "_gp.csv").write(gp_sol);

    return 0;
}