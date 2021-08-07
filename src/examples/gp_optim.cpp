#include <iostream>

#include <Eigen/Core>
#include <gp_manifold/GaussianProcess.hpp>
#include <gp_manifold/optimization/IpoptOptimizer.hpp>
#include <utils_cpp/UtilsCpp.hpp>

using namespace gp_manifold;
using namespace kernel_lib;

struct ParamsExp {
    struct kernel : public defaults::kernel {
        PARAM_SCALAR(double, sf, -0.2);
        PARAM_SCALAR(double, sn, -0.4);
    };

    struct exp_sq : public defaults::exp_sq {
        PARAM_SCALAR(double, l, -0.3566);
    };

    struct exp_sq_full : public defaults::exp_sq_full {
    };

    struct gaussian : public defaults::gaussian {
    };
};

int main(int argc, char** argv)
{
    // data
    int num_sample = 50, res = 200;
    double range = 10;

    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(res, -range, range),
                    y = x.array().sin(),
                    sample_x = range * Eigen::VectorXd::Random(num_sample),
                    sample_y = sample_x.array().sin() + 0.2 * Eigen::VectorXd::Random(num_sample).array();

    // Gaussian Process
    using Kernel_t = kernels::SquaredExp<ParamsExp>;
    using GP_t = GaussianProcess<ParamsExp, Kernel_t, SmartPtr<optimization::IpoptOptimizer>>;
    GP_t gp;

    // Set training point and target
    gp.setSamples(sample_x).setTarget(sample_y).update();

    // Save GP solution
    utils_cpp::FileManager io_manager;
    io_manager.setFile("rsc/solutions/temp_gp.csv").write("X", x, "Y", y, "sample_x", sample_x, "sample_y", sample_y, "GP", gp.multiEval(x));

    // gp.check();

    // Optimize gp
    if (gp.optimize())
        std::cout << "Model OPTIMIZED" << std::endl;
    else
        std::cout << "Model NOT optimized" << std::endl;

    io_manager.append("OPT", gp.multiEval(x));

    return 0;
}