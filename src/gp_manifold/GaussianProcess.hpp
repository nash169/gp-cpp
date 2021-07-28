#ifndef GP_MANIFOLD_GAUSSIAN_PROCESS
#define GP_MANIFOLD_GAUSSIAN_PROCESS

#include <kernel_lib/Kernel.hpp>

namespace gp_manifold {
    using namespace kernel_lib;

    template <typename Params, typename Kernel>
    class GaussianProcess : public utils::Expansion<Params, Kernel> {
    public:
        GaussianProcess()
        {
        }

        GaussianProcess& setReference(const Eigen::MatrixXd& reference) override
        {
            // utils::Expansion<Params, Kernel>::setReference(reference);

            // return *this;
            return static_cast<GaussianProcess&>(utils::Expansion<Params, Kernel>::setReference(reference));
        }

        GaussianProcess& setTarget(const Eigen::VectorXd& target)
        {
            _target = target;

            return *this;
        }

        GaussianProcess& update()
        {
            Eigen::MatrixXd G = utils::Expansion<Params, Kernel>::_kernel
                                    .gram(utils::Expansion<Params, Kernel>::_reference, utils::Expansion<Params, Kernel>::_reference);

            Eigen::LLT<Eigen::MatrixXd, Eigen::Upper> U = G.selfadjointView<Eigen::Upper>().llt();

            // std::cout << U.solve(_target) << std::endl;

            utils::Expansion<Params, Kernel>::setParams(U.solve(_target));

            return *this;
        }

    protected:
        // Training points
        Eigen::VectorXd _target;

        // Gaussian distribution
        utils::Gaussian<Params, kernels::SquaredExpFull<Params>> _gauss;

        double likelihood()
        {
            return _gauss.log(_target);
        }

        auto likelihoodGrad()
        {
            Eigen::VectorXd grad(_kernel.sizeParams()), gradKernel = _kernel.gradParams(_target, _target);
            Eigen::MatrixXd gradGauss = Eigen::Map<Eigen::MatrixXd>(_gauss.gradParams(_target).segment(0, _target.size() * _target.size()).data(), _target.size(), _target.size());

            for (size_t i = 0; i < _kernel.sizeParams(); i++)
                grad(i) = 0.5 * (gradGauss * gradKernel(i)).trace();

            return grad;
        }
    };
} // namespace gp_manifold

#endif // GP_MANIFOLD_GAUSSIAN_PROCESS