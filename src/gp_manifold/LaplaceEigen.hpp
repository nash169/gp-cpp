#ifndef GP_MANIFOLD_LAPLACE_EIGEN
#define GP_MANIFOLD_LAPLACE_EIGEN

#include <kernel_lib/Kernel.hpp>

namespace gp_manifold {
    using namespace kernel_lib;

    template <typename Params, typename Kernel>
    class LaplaceEigen : public utils::EigenFunction<utils::Expansion<Params, Kernel>> {
    public:
        LaplaceEigen()
        {
        }

        // Override operator()
        inline double operator()(const Eigen::Matrix<double, Eigen::Dynamic, 1>& x, const double& eigenvalue) const override
        {
            return utils::EigenFunction<utils::Expansion<Params, Kernel>>::_eigen_pair.at(eigenvalue)(x);
        }

    protected:
    };
} // namespace gp_manifold

#endif // GP_MANIFOLD_LAPLACE_EIGEN