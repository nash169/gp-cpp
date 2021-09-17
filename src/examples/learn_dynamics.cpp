#include <iostream>

#include <kernel_lib/Kernel.hpp>

#include <utils_cpp/UtilsCpp.hpp>

#include <gp_manifold/SlepcSolver.hpp>

using namespace kernel_lib;
using namespace gp_manifold;

struct Params {
    struct kernel : public defaults::kernel {
    };

    struct exp_sq : public defaults::exp_sq {
        PARAM_SCALAR(double, l, -2.30259);
    };
};

struct Embedding {
    Embedding()
    {
        Eigen::Matrix<double, 1, 2> x;

        // Expansion 1
        x << -25, 25;
        psi1.setSamples(x);
        psi1.setWeights(tools::makeVector(1));
        psi1.kernel().setParams(Eigen::Vector3d(0, -1e-2, 2));

        // Expansion 2
        x << 25, 25;
        psi2.setSamples(x);
        psi2.setWeights(tools::makeVector(1));
        psi2.kernel().setParams(Eigen::Vector3d(0, -1e-2, 1.5));

        // Expansion 3
        x << 25, -25;
        psi3.setSamples(x);
        psi3.setWeights(tools::makeVector(1));
        psi3.kernel().setParams(Eigen::Vector3d(0, -1e-2, 0.5));
    }

    Eigen::Vector3d operator()(const Eigen::Vector2d& x)
    {
        return Eigen::Vector3d(x(0), x(1), psi1(x) + psi2(x) + psi3(x));
    }

    Eigen::Matrix<double, Eigen::Dynamic, 3> multiEval(const Eigen::Matrix<double, Eigen::Dynamic, 2>& x)
    {
        Eigen::Matrix<double, Eigen::Dynamic, 3> embedding(x.rows(), 3);
        embedding.col(0) = x.col(0);
        embedding.col(1) = x.col(1);
        embedding.col(2) = psi1.multiEval(x) + psi2.multiEval(x) + psi3.multiEval(x);

        return embedding;
    }

    Eigen::Matrix<double, 3, 2> grad(const Eigen::Vector2d& x)
    {
        Eigen::Matrix<double, 3, 2> g;
        g.row(0) << 1, 0;
        g.row(1) << 0, 1;
        g.row(2) << (psi1.grad(x) + psi2.grad(x) + psi3.grad(x)).transpose();

        return g;
    }

    Eigen::Matrix<double, Eigen::Dynamic, 6> multiGrad(const Eigen::Matrix<double, Eigen::Dynamic, 2>& x)
    {
        Eigen::Matrix<double, Eigen::Dynamic, 6> g(x.rows(), 6);
        g.col(0).setConstant(1);
        g.col(1).setConstant(0);
        g.col(2).setConstant(0);
        g.col(3).setConstant(1);
        g.block(0, 4, g.rows(), 2) = psi1.multiGrad(x) + psi2.multiGrad(x) + psi3.multiGrad(x);

        return g;
    }

    using Expansion_t = utils::Expansion<Params, kernels::SquaredExp<Params>>;
    Expansion_t psi1, psi2, psi3;
};

struct Dynamics {
    Dynamics()
    {
        A = -Eigen::Matrix2d::Identity();
    }

    Eigen::Vector2d operator()(const Eigen::Vector2d& x)
    {
        return A * x;
    }

    Eigen::Matrix<double, 2, 2> A;
};

struct Potential {
    Potential()
    {
        P = Eigen::Matrix3d::Identity();
    }

    double operator()(const Eigen::Vector3d& x)
    {
        return x.transpose() * P * x;
    }

    Eigen::Matrix3d P;
};

int main(int argc, char** argv)
{
    double box[] = {-50, 50, -50, 50};
    size_t resolution = 100, num_samples = resolution * resolution, dim = sizeof(box) / sizeof(*box) / 2;

    // Data
    Eigen::MatrixXd X = Eigen::RowVectorXd::LinSpaced(resolution, box[0], box[1]).replicate(resolution, 1),
                    Y = Eigen::VectorXd::LinSpaced(resolution, box[2], box[3]).replicate(1, resolution),
                    x_test(num_samples, dim);
    x_test << Eigen::Map<Eigen::VectorXd>(X.data(), X.size()), Eigen::Map<Eigen::VectorXd>(Y.data(), Y.size());

    // Dynamics
    Eigen::Matrix<double, Eigen::Dynamic, 3> field3d(x_test.rows(), 3);
    for (size_t i = 0; i < field3d.rows(); i++) {
        Eigen::Vector2d x = x_test.row(i);
        field3d.row(i) = Embedding().grad(x) * Dynamics()(x);
    }

    // Samples
    Eigen::Matrix<double, Eigen::Dynamic, 3> samples = Embedding().multiEval(x_test);

    // Geometric Laplacian
    double eps = 2 * std::pow(std::exp(-2.30259), 2);
    size_t nn = 50, skip_diag = 1;
    utils::Graph graph;

    Eigen::SparseMatrix<double, Eigen::RowMajor> L = graph.kNearestWeighted(samples, kernels::SquaredExp<Params>(), nn, skip_diag);
    std::vector<Eigen::Triplet<double>> tripletList;
    for (size_t i = 0; i < L.rows(); i++)
        tripletList.push_back(Eigen::Triplet<double>(i, i, 1 / L.row(i).sum()));
    Eigen::SparseMatrix<double, Eigen::RowMajor> D(num_samples, num_samples);
    D.setFromTriplets(tripletList.begin(), tripletList.end());
    L = D * L * D;
    tripletList.clear();
    for (size_t i = 0; i < L.rows(); i++)
        tripletList.push_back(Eigen::Triplet<double>(i, i, 1 / L.row(i).sum()));
    D.setFromTriplets(tripletList.begin(), tripletList.end());
    L = (Eigen::MatrixXd::Identity(num_samples, num_samples) - D * L) / eps / 4;

    // Create Slepc solver
    int nev = 100;
    SlepcSolver solver(argc, argv);

    solver.operators(std::make_unique<PetscMatrix>(L)) // Set opeartors
        .target(0.0) // target eigenvalue
        .spectrum(EPS_TARGET_REAL) // Select smallest eigenvalues
        .modes(nev) // Number of requested eigenvalues
        .spectralTransformation(STSINVERT)
        .solve(); // Solve the problem

    // Number of eigenvalues converged
    std::cout << "Number of converged eigenvalues: " << solver.converged() << std::endl;

    // Save eigenvalues & eigenvector
    Eigen::VectorXd eigs(nev);
    Eigen::MatrixXd vecs(num_samples, nev);

    for (size_t i = 0; i < nev; i++) {
        eigs(i) = solver.eigenvalue(i);
        vecs.col(i) = solver.eigenvector(i);
    }

    // Potential
    Eigen::VectorXd fun(samples.rows());
    for (size_t i = 0; i < samples.rows(); i++) {
        Eigen::Vector3d x = samples.row(i);
        Potential()(x);
    }

    // Solution
    utils_cpp::FileManager io_manager;
    io_manager.setFile("rsc/solutions/kernel_eval.csv");
    io_manager.write("SAMPLES", samples, "FUNCTION", fun, "EIGVEC", vecs);

    return 0;
}