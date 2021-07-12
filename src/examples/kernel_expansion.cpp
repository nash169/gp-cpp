#include <iostream>
#include <random>
#include <sstream>

#include <Eigen/Core>

#include <magnum_dynamics/MagnumApp.hpp>
#include <utils_cpp/UtilsCpp.hpp>

#include <gp_manifold/GaussianProcess.hpp>
#include <gp_manifold/LaplaceEigen.hpp>

using namespace gp_manifold;
using namespace kernel_lib;

struct ParamsExp {
    struct kernel : public defaults::kernel {
        PARAM_SCALAR(double, sf, 0);
        PARAM_SCALAR(double, sn, -5);
    };

    struct exp_sq : public defaults::exp_sq {
        PARAM_SCALAR(double, l, std::log(0.7));
    };
};

struct ParamsRiemann {
    struct kernel : public defaults::kernel {
        PARAM_SCALAR(double, sf, 0);
        PARAM_SCALAR(double, sn, -5);
    };

    struct riemann_exp_sq : public defaults::riemann_exp_sq {
        PARAM_SCALAR(double, l, std::log(0.1));
    };
};

constexpr int MIN = 0;
constexpr int MAX = 1595;

constexpr int RAND_NUMS_TO_GENERATE = 1000;

int main(int argc, char** argv)
{
    // Load nodes, eigenvectors and eigenvalues
    utils_cpp::FileManager io_manager;
    Eigen::MatrixXd vertices = io_manager.setFile("rsc/sol/sphere.mesh").read<Eigen::MatrixXd>("vertices", 3);

    size_t num_eig = 30;
    Eigen::MatrixXd eigenvectors(vertices.rows(), num_eig);

    for (size_t i = 0; i < num_eig; i++) {
        std::stringstream file_path;
        file_path << "rsc/sol/sphere_mode_" << i;
        eigenvectors.col(i) = io_manager.setFile(file_path.str()).read<Eigen::MatrixXd>("", 5);
    }

    Eigen::VectorXd eigenvalues = io_manager.setFile("rsc/sol/sphere_eigs").read<Eigen::MatrixXd>();

    // Extract reference and ground truth
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_int_distribution<int> distr(MIN, MAX);

    Eigen::MatrixXd nodes = io_manager.setFile("rsc/data/mesh_vertices.csv").read<Eigen::MatrixXd>(),
                    indices = io_manager.setFile("rsc/data/mesh_faces.csv").read<Eigen::MatrixXd>(),
                    reference(RAND_NUMS_TO_GENERATE, nodes.cols());

    Eigen::VectorXd ground_truth = io_manager.setFile("rsc/data/ground_truth.csv").read<Eigen::MatrixXd>(),
                    target(RAND_NUMS_TO_GENERATE);

    for (size_t i = 0; i < RAND_NUMS_TO_GENERATE; i++) {
        int index = distr(eng);
        reference.row(i) = nodes.row(index);
        target(i) = ground_truth(index);
    }

    using Kernel_t = kernels::SquaredExp<ParamsExp>;
    using Expansion_t = utils::Expansion<ParamsExp, Kernel_t>;

    // Gaussian Process
    using GP_t = GaussianProcess<ParamsExp, Kernel_t>;
    GP_t gp;

    gp.setReference(reference).setTarget(target).update();

    Eigen::VectorXd gp_sol = gp.multiEval(nodes);

    io_manager.setFile("rsc/data/gp_sol.csv").write(gp_sol);

    // Riemann Gaussian Process
    using Riemann_t = kernels::RiemannSqExp<ParamsRiemann, Expansion_t>;
    using RGP_t = GaussianProcess<ParamsRiemann, Riemann_t>;
    RGP_t rgp;

    for (size_t i = 1; i < num_eig; i++) {
        // Create eigenfunction
        Expansion_t f;

        // Set manifold sampled points and weights
        f.setReference(vertices).setParams(eigenvectors.col(i));

        // Add eigen-pair to Riemann kernel
        rgp.kernel().addPair(eigenvalues(i), f);
    }

    rgp.setReference(reference).setTarget(target).update();

    Eigen::VectorXd rgp_sol = rgp.multiEval2(nodes);

    io_manager.setFile("rsc/data/rgp_sol.csv").write(rgp_sol);

    return 0;
}