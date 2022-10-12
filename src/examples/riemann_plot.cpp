#include <iostream>
#include <random>
#include <sstream>

#include <gp_manifold/GaussianProcess.hpp>
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
    std::string mesh_name = "torus",
                mesh_ext = "mesh";
    FileManager io_manager;

    // Load mesh nodes and indices
    Eigen::MatrixXd vertices = io_manager.setFile("outputs/modes/fem_" + mesh_name + "_mesh.000000").read<Eigen::MatrixXd>("vertices", 3);

    // Load eigenvalues & eigenvectors
    Eigen::VectorXd eigenvalues = io_manager.setFile("outputs/modes/fem_" + mesh_name + "_eigs.000000").read<Eigen::MatrixXd>("eigs", 2);
    Eigen::MatrixXd eigenvectors = io_manager.setFile("outputs/modes/fem_" + mesh_name + "_modes.000000").read<Eigen::MatrixXd>("modes", 2);

    // Riemannian Gaussian Process
    using Kernel_t = kernels::SquaredExp<ParamsExp>;
    using Expansion_t = utils::Expansion<ParamsExp, Kernel_t>;

    using Riemann_t = kernels::RiemannSqExp<ParamsRiemann, Expansion_t>;
    using RiemannExpansion_t = utils::Expansion<ParamsRiemann, Riemann_t>;

    RiemannExpansion_t phi;

    // Set kernel eigen pairs
    int num_modes = (argc > 2) ? std::stoi(argv[2]) : 10;

    for (size_t i = 1; i < num_modes; i++) {
        // Create eigenfunction
        Expansion_t f;

        // Set manifold sampled points and weights
        f.setSamples(vertices).setWeights(eigenvectors.col(i));

        // Add eigen-pair to Riemann kernel
        phi.kernel().addPair(eigenvalues(i), f);
    }

    // Test space
    double box[] = {0, 2 * M_PI, 0, 2 * M_PI};
    size_t resolution = 100, num_samples = resolution * resolution, dim = sizeof(box) / sizeof(*box) / 2;

    // Data
    Eigen::MatrixXd X = Eigen::RowVectorXd::LinSpaced(resolution, box[0], box[1]).replicate(resolution, 1),
                    Y = Eigen::VectorXd::LinSpaced(resolution, box[2], box[3]).replicate(1, resolution),
                    x_chart(num_samples, dim), x_ref(1, dim + 1);

    // Embedding
    auto sphere_embed = [](const Eigen::MatrixXd& x) {
        Eigen::MatrixXd e(x.rows(), 3);
        e.col(0) = x.col(0).array().sin() * x.col(1).array().cos();
        e.col(1) = x.col(0).array().sin() * x.col(1).array().sin();
        e.col(2) = x.col(0).array().cos();
        return e;
    };

    auto torus_embed = [](const Eigen::MatrixXd& x) {
        double a = 1, b = 3;
        Eigen::MatrixXd e(x.rows(), 3);
        e.col(0) = (b + a * x.col(0).array().cos()) * x.col(1).array().cos();
        e.col(1) = (b + a * x.col(0).array().cos()) * x.col(1).array().sin();
        e.col(2) = a * x.col(0).array().sin();
        return e;
    };

    // Random number generator
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_int_distribution<int> distr(0, x_chart.rows());

    size_t ref_id = 7522; // distr(eng);
    std::cout << ref_id << std::endl;
    x_chart << Eigen::Map<Eigen::VectorXd>(X.data(), X.size()), Eigen::Map<Eigen::VectorXd>(Y.data(), Y.size());
    x_ref << torus_embed(x_chart.row(ref_id));

    // Set training point and target
    phi.setSamples(x_ref).setWeights(tools::makeVector(1));

    // Save GP solution
    io_manager
        .setFile("outputs/chart_riemann.csv")
        .write(x_chart);
    io_manager
        .setFile("outputs/embed_riemann.csv")
        .write(torus_embed(x_chart));
    io_manager
        .setFile("outputs/eval_riemann.csv")
        .write(phi.multiEval2(torus_embed(x_chart)));

    return 0;
}