#ifndef GP_MANIFOLD_GAUSSIAN_PROCESS
#define GP_MANIFOLD_GAUSSIAN_PROCESS

#include <kernel_lib/Kernel.hpp>

#include "gp_manifold/optimization/AbstractOptimizer.hpp"

namespace gp_manifold {
    using namespace kernel_lib;

    template <typename Params, typename Kernel, typename Optimizer = optimization::AbstractOptimizer>
    class GaussianProcess : public utils::Expansion<Params, Kernel> {
    public:
        GaussianProcess() : _opt()
        {
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
            params_gauss.segment(0, _y.size() * _y.size()) = Eigen::Map<Eigen::VectorXd>(_k.gram(_y, _y).data(), _y.size() * _y.size());

            // Set Gaussian params
            _g.setParams(params_gauss);
        }

        double likelihood(const Eigen::VectorXd& params)
        {
            // Set params likelihood (no mean optimization at the moment)
            setLikelihood(params);

            return _g.log(_y);
        }

        auto likelihoodGrad(const Eigen::VectorXd& params)
        {
            // Set params likelihood (no mean optimization at the moment)
            setLikelihood(params);

            // Params gradient of the kernel
            Eigen::VectorXd grad(_k.sizeParams()), gradKernel = _k.gradParams(_y, _y);

            // Params gradient of the Gaussian distribution
            Eigen::MatrixXd gradGauss = Eigen::Map<Eigen::MatrixXd>(_g.gradParams(_y).segment(0, _y.size() * _y.size()).data(), _y.size(), _y.size());

            // Total gradient
            for (size_t i = 0; i < _k.sizeParams(); i++)
                grad(i) = 0.5 * (gradGauss * gradKernel(i)).trace();

            return grad;
        }

        void optimize()
        {
            // For passing methods of a class: https://stackoverflow.com/questions/7582546/using-generic-stdfunction-objects-with-member-functions-in-one-class
            // https://stackoverflow.com/questions/26331628/reference-to-non-static-member-function-must-be-called

            // double (GaussianProcess<Params, Kernel, Optimizer>::*function)(const Eigen::VectorXd&) = &GaussianProcess<Params, Kernel, Optimizer>::likelihood;
            _opt.setObjective(std::bind(&GaussianProcess<Params, Kernel, Optimizer>::likelihood, this, std::placeholders::_1));
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