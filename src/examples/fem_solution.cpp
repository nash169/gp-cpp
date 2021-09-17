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
    std::string mesh_name = "armadillo", mesh_ext = "mesh";
    utils_cpp::FileManager io_manager;

    // Load mesh nodes and indices
    Eigen::MatrixXd vertices = io_manager.setFile("rsc/modes/fem_" + mesh_name + "_mesh.000000").read<Eigen::MatrixXd>("vertices", 3),
                    indices = io_manager.read<Eigen::MatrixXd>("elements", 2);

    // Load eigenvalues & eigenvectors
    Eigen::VectorXd eigenvalues = io_manager.setFile("rsc/modes/fem_" + mesh_name + "_eigs.000000").read<Eigen::MatrixXd>("eigs", 2);
    Eigen::MatrixXd eigenvectors = io_manager.setFile("rsc/modes/fem_" + mesh_name + "_modes.000000").read<Eigen::MatrixXd>("modes", 2);

    // Load ground truth, target and relative nodes
    Eigen::MatrixXd nodes = io_manager.setFile("rsc/truth/" + mesh_name + "_vertices.csv").read<Eigen::MatrixXd>(),
                    faces = io_manager.setFile("rsc/truth/" + mesh_name + "_faces.csv").read<Eigen::MatrixXd>(),
                    reference = io_manager.setFile("rsc/truth/" + mesh_name + "_reference.csv").read<Eigen::MatrixXd>();

    Eigen::VectorXd ground_truth = io_manager.setFile("rsc/truth/" + mesh_name + "_truth.csv").read<Eigen::MatrixXd>(),
                    target = io_manager.setFile("rsc/truth/" + mesh_name + "_target.csv").read<Eigen::MatrixXd>();

    // Riemannian Gaussian Process
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
    io_manager.setFile("rsc/solutions/fem_" + mesh_name + "_rgp.csv").write(rgp_sol);

    return 0;
}