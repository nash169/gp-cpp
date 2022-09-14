#include <Eigen/Core>
#include <iostream>
#include <kernel_lib/kernels/SquaredExp.hpp>
#include <random>
#include <utils_lib/FileManager.hpp>
#include <vector>

using namespace kernel_lib;
using namespace utils_lib;

struct Params {
    struct kernel : public defaults::kernel {
    };

    struct exp_sq : public defaults::exp_sq {
    };
};

Eigen::MatrixXd uniformSampleSphere(const Eigen::VectorXd& center, const double& radius = 1, const size_t& num_points = 1000)
{
    Eigen::MatrixXd X = Eigen::MatrixXd::Random(num_points, center.size());
    X.rowwise().normalize();
    X.rowwise() += center.transpose();

    return X;
}

int main(int argc, char const* argv[])
{
    kernels::SquaredExp<Params> k;

    // Sphere uniform sampling
    Eigen::MatrixXd X = uniformSampleSphere(Eigen::VectorXd::Zero(3));

    // Random number generator
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_int_distribution<int> distr(0, X.rows());

    // Kernel evaluation
    Eigen::MatrixXd Xi = X.row(distr(eng));
    Eigen::VectorXd f = k.gram(Xi, X).transpose();

    // Save
    FileManager man("rsc/tangent.csv");
    man.write("reference", Xi, "samples", X, "kernel", f);

    // Plot
    // std::vector<double> vec1(X.col(0).data(), X.col(0).data() + X.rows()),
    //     vec2(X.col(1).data(), X.col(1).data() + X.rows()),
    //     vec3(X.col(2).data(), X.col(2).data() + X.rows());
    // plt::scatter(vec1, vec2, vec3);
    // plt::show();

    return 0;
}
