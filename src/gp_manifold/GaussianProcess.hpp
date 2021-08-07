#ifndef GP_MANIFOLD_GAUSSIAN_PROCESS
#define GP_MANIFOLD_GAUSSIAN_PROCESS

#include <memory>
#include <type_traits>

#include <kernel_lib/Kernel.hpp>
#include <utils_cpp/UtilsCpp.hpp>

#include "gp_manifold/optimization/AbstractOptimizer.hpp"

namespace gp_manifold {
    using namespace kernel_lib;

    template <typename Params, typename Kernel, typename Optimizer = std::unique_ptr<optimization::AbstractOptimizer>>
    class GaussianProcess : public utils::Expansion<Params, Kernel> {
    public:
        GaussianProcess()
        {
            _opt = Optimizer(new typename std::remove_reference<decltype(*_opt)>::type);
        }

        GaussianProcess& setSamples(const Eigen::MatrixXd& x) override
        {
            // utils::Expansion<Params, Kernel>::setReference(reference);

            // return *this;
            return static_cast<GaussianProcess&>(utils::Expansion<Params, Kernel>::setSamples(x));
        }

        GaussianProcess& setTarget(const Eigen::VectorXd& y)
        {
            _y = y;

            return *this;
        }

        GaussianProcess& update()
        {
            Eigen::MatrixXd G = _k.gram(_x, _x);

            Eigen::LLT<Eigen::MatrixXd, Eigen::Upper> U = G.selfadjointView<Eigen::Upper>().llt();

            // std::cout << U.solve(_y) << std::endl;

            setWeights(U.solve(_y));

            return *this;
        }

        void setLikelihood(const Eigen::VectorXd& params)
        {
            // Set covariance function parameters (no mean optimization at the moment)
            _k.setParams(params); // This is setting signal and noise variance as well

            // Set zero mean and covariance matrix for the Gaussian distribution
            Eigen::VectorXd params_gauss = Eigen::VectorXd::Zero(_y.size() * _y.size() + _y.size());

            // Set gram matrix
            params_gauss.segment(0, _y.size() * _y.size()) = Eigen::Map<Eigen::VectorXd>(_k.gram(_x, _x).data(), _y.size() * _y.size());

            // Set Gaussian params
            _g.setParams(params_gauss);
        }

        double likelihood(const Eigen::VectorXd& params)
        {
            // Set params likelihood (no mean optimization at the moment)
            setLikelihood(params);

            return -_g.log(_y);
        }

        auto likelihoodGrad(const Eigen::VectorXd& params)
        {
            // Set params likelihood (no mean optimization at the moment)
            setLikelihood(params);

            // Params gradient of the kernel
            Eigen::MatrixXd gradKernel = _k.gramGradParams(_x, _x);

            // Params gradient of the Gaussian distribution (only covariance)
            Eigen::MatrixXd gradGauss = Eigen::Map<Eigen::MatrixXd>(_g.gradParams(_y).segment(0, _y.size() * _y.size()).data(), _y.size(), _y.size());

            // Total gradient
            Eigen::VectorXd grad(_k.sizeParams());
            for (size_t i = 0; i < _k.sizeParams(); i++)
                grad(i) = (gradGauss * Eigen::Map<Eigen::MatrixXd>(gradKernel.col(i).data(), _y.size(), _y.size())).trace();

            return -grad;
        }

        bool optimize()
        {
            // For passing methods of a class: https://stackoverflow.com/questions/7582546/using-generic-stdfunction-objects-with-member-functions-in-one-class
            // https://stackoverflow.com/questions/26331628/reference-to-non-static-member-function-must-be-called

            // double (GaussianProcess<Params, Kernel, Optimizer>::*function)(const Eigen::VectorXd&) = &GaussianProcess<Params, Kernel, Optimizer>::likelihood;
            (*_opt)
                .setDimension(this->sizeParams())
                .setStartingPoint(this->params())
                .setObjective(std::bind(&GaussianProcess<Params, Kernel, Optimizer>::likelihood, this, std::placeholders::_1))
                .setObjectiveGradient(std::bind(&GaussianProcess<Params, Kernel, Optimizer>::likelihoodGrad, this, std::placeholders::_1));

            return _opt->optimize();
        }

        bool check()
        {
            std::cout << "Number of params: " << this->sizeParams() << std::endl;

            utils_cpp::DerivativeChecker checker(this->sizeParams());

            return checker.checkGradient(std::bind(&GaussianProcess<Params, Kernel, Optimizer>::likelihood, this, std::placeholders::_1),
                std::bind(&GaussianProcess<Params, Kernel, Optimizer>::likelihoodGrad, this, std::placeholders::_1));
        }

        // double sigma(const Eigen::VectorXd& x)
        // {
        //     Eigen::MatrixXd G = _k.gram(_x, _x);

        //     Eigen::LLT<Eigen::MatrixXd, Eigen::Upper> U = G.selfadjointView<Eigen::Upper>().llt();

        //     Eigen::VectorXd y = _k(_x, x);

        //     return _k(x, x) - y.transpose() * U.solve(y);
        // }

        // Eigen::VectorXd multiSigma(const Eigen::MatrixXd& x)
        // {
        //     Eigen::MatrixXd G = _k.gram(_x, _x);

        //     Eigen::LLT<Eigen::MatrixXd, Eigen::Upper> U = G.selfadjointView<Eigen::Upper>().llt();

        //     Eigen::VectorXd sigma(x.rows()), y = _k(_x, x);

        //     for (size_t i = 0; i < x.rows(); i++)
        //     sigma(i) = _k(x.row(i), x) - y.transpose() * U.solve(y);

        //     return sigma;
        // }

        Eigen::VectorXd params() const
        {
            return _k.params();
        }

        GaussianProcess& setParams(const Eigen::VectorXd& params)
        {
            _k.setParams(params);
        }

        size_t sizeParams() const
        {
            return _k.sizeParams();
        }

    protected:
        // Lift samples to GaussianProcess
        using utils::Expansion<Params, Kernel>::_x;

        // Lift kernel to GaussianProcess
        using utils::Expansion<Params, Kernel>::_k;

        // Make protected setWeights method
        using utils::Expansion<Params, Kernel>::setWeights;

        // Training points
        Eigen::VectorXd _y;

        // Gaussian distribution
        utils::Gaussian<Params, kernels::SquaredExpFull<Params>> _g;

        // Optimizer
        Optimizer _opt;
    };
} // namespace gp_manifold

#endif // GP_MANIFOLD_GAUSSIAN_PROCESS