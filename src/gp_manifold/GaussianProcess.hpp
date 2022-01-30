#ifndef GP_MANIFOLD_GAUSSIAN_PROCESS
#define GP_MANIFOLD_GAUSSIAN_PROCESS

#include <memory>
#include <type_traits>

#include <kernel_lib/Kernel.hpp>
#include <utils_lib/DerivativeChecker.hpp>

#include "gp_manifold/optimization/AbstractOptimizer.hpp"

using namespace kernel_lib;
using namespace utils_lib;

namespace gp_manifold {

    /* Defaults parameters for the normal distribution (those cannot be initialize in the main) */
    struct GaussianParams {
        struct gaussian : public defaults::gaussian {
        };

        struct kernel : public defaults::kernel {
        };

        struct exp_sq_full : public defaults::exp_sq_full {
        };
    };

    // namespace defaults {
    //     struct gaussian_process {
    //     };
    // } // namespace defaults

    template <typename Params, typename Kernel, typename Optimizer = std::unique_ptr<optimization::AbstractOptimizer>>
    class GaussianProcess : public utils::Expansion<Params, Kernel> {
    public:
        GaussianProcess()
        {
            _opt = Optimizer(new typename std::remove_reference<decltype(*_opt)>::type);
        }

        /* Set samples */
        GaussianProcess& setSamples(const Eigen::MatrixXd& x) override
        {
            // utils::Expansion<Params, Kernel>::setReference(reference);

            // return *this;
            return static_cast<GaussianProcess&>(utils::Expansion<Params, Kernel>::setSamples(x));
        }

        /* Set target */
        GaussianProcess& setTarget(const Eigen::VectorXd& y)
        {
            _y = y;

            return *this;
        }

        /* Update gaussian process model (this automatically perform the cholesky decomposition within the Gaussian) */
        GaussianProcess& update()
        {
            // Set Gaussian covariance to automatically perform the cholesky decomposition (within SqExpFull)
            // for now we use setGaussian because the mean is set to zero (fix this)
            // put some condition to check if _x is present
            setGaussian();

            // Eigen::MatrixXd G = _k.gram(_x, _x);

            // Eigen::LLT<Eigen::MatrixXd, Eigen::Upper> U = G.selfadjointView<Eigen::Upper>().llt();

            // setWeights(U.solve(_y));

            setWeights(_g.template covariance<Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>, Eigen::Lower>>().solve(_y));

            return *this;
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

            DerivativeChecker checker(this->sizeParams());

            return checker.checkGradient(std::bind(&GaussianProcess<Params, Kernel, Optimizer>::likelihood, this, std::placeholders::_1),
                std::bind(&GaussianProcess<Params, Kernel, Optimizer>::likelihoodGrad, this, std::placeholders::_1));
        }

        double sigma(const Eigen::VectorXd& x)
        {
            // Eigen::MatrixXd G = _k.gram(_x, _x);

            // Eigen::LLT<Eigen::MatrixXd, Eigen::Upper> U = G.selfadjointView<Eigen::Upper>().llt();

            Eigen::VectorXd y = _k.gram(_x, Eigen::MatrixXd(x.transpose()));

            return _k(x, x) - y.transpose() * _g.template covariance<Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>, Eigen::Lower>>().solve(y);
        }

        Eigen::VectorXd multiSigma(const Eigen::MatrixXd& x)
        {
            // Eigen::MatrixXd G = _k.gram(_x, _x);

            // Eigen::LLT<Eigen::MatrixXd, Eigen::Upper> U = G.selfadjointView<Eigen::Upper>().llt();

            Eigen::VectorXd sigma(x.rows());
            Eigen::MatrixXd y = _k.gram(_x, x);

#pragma omp parallel for
            for (size_t i = 0; i < x.rows(); i++)
                sigma(i) = _k(x.row(i), x.row(i)) - y.col(i).transpose() * _g.template covariance<Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>, Eigen::Lower>>().solve(y.col(i));

            return sigma;
        }

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
        void setGaussian()
        {
            // Set zero mean and covariance matrix for the Gaussian distribution
            Eigen::VectorXd params_gauss = Eigen::VectorXd::Zero(_y.size() * _y.size() + _y.size());

            // Set gram matrix
            params_gauss.segment(0, _x.rows() * _x.rows()) = Eigen::Map<Eigen::VectorXd>(_k.gram(_x, _x).data(), _x.rows() * _x.rows());

            // Set Gaussian params
            _g.setParams(params_gauss);
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

        // Lift samples to GaussianProcess
        using utils::Expansion<Params, Kernel>::_x;

        // Lift kernel to GaussianProcess
        using utils::Expansion<Params, Kernel>::_k;

        // Lift and make protected setWeights method
        using utils::Expansion<Params, Kernel>::setWeights;
        // Logic behind using: https://stackoverflow.com/questions/20790932/what-is-the-logic-behind-the-using-keyword-in-c

        // Training points
        Eigen::VectorXd _y;

        // Gaussian distribution
        utils::Gaussian<GaussianParams, kernels::SquaredExpFull<GaussianParams>> _g;

        // Optimizer
        Optimizer _opt;
    };
} // namespace gp_manifold

#endif // GP_MANIFOLD_GAUSSIAN_PROCESS