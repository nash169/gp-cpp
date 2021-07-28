#ifndef GP_MANIFOLD_ABSTRACT_OPTIMIZER
#define GP_MANIFOLD_ABSTRACT_OPTIMIZER

namespace gp_manifold {
    namespace optimization {
        class AbstractOptimizer {
        public:
            AbstractOptimizer()
            {
            }

            ~AbstractOptimizer()
            {
            }

            template <typename Function>
            AbstractOptimizer& setFunction(const Function& f)
            {
            }

        protected:
        };
    } // namespace optimization
} // namespace gp_manifold

#endif // GP_MANIFOLD_ABSTRACT_OPTIMIZER