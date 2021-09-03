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

    // Loads eigenvalues and eigenvectors
    Eigen::VectorXd eigenvalues = io_manager.setFile("rsc/modes/diffusion_" + mesh_name + "_eigs.000000").read<Eigen::MatrixXd>("eigs", 2);
    Eigen::MatrixXd eigenvectors = io_manager.setFile("rsc/modes/diffusion_" + mesh_name + "_modes.000000").read<Eigen::MatrixXd>("modes", 2);

    // Riemann Gaussian Process
    using Kernel_t = kernels::SquaredExp<ParamsExp>;
    using Expansion_t = utils::Expansion<ParamsExp, Kernel_t>;
    using Riemann_t = kernels::RiemannSqExp<ParamsRiemann, Expansion_t>;
    using RGP_t = GaussianProcess<ParamsRiemann, Riemann_t>;
    RGP_t rgp;

    // Set kernel eigen pairs
    int num_modes = 100;

    for (size_t i = 1; i < num_modes; i++) {
        // Create eigenfunction
        Expansion_t f;

        // Set manifold sampled points and weights
        f.setSamples(vertices).setWeights(eigenvectors.col(i));

        // Add eigen-pair to Riemann kernel
        rgp.kernel().addPair(eigenvalues(i), f);
    }

    // Set training point and target
    rgp.setSamples(reference).setTarget(target).update();

    // Evaluation of the GP on all the mesh points
    Eigen::VectorXd rgp_sol = rgp.multiEval2(nodes);

    // Save GP solution
    io_manager.setFile("rsc/solutions/diffusion_" + mesh_name + "_rgp.csv").write(rgp_sol);

    return 0;
}