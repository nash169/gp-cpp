#ifndef GPMANIFOLD_EIGEN_HPP
#define GPMANIFOLD_EIGEN_HPP

#include <Eigen/Core>
#include <mfem.hpp> // mfem/linalg/operator.hpp

namespace gp_manifold {
    class EigenEigenSolver {
    public:
        /**
         * @brief Initializes the FrankaHW class to be fully operational. This involves parsing required
         * configurations from the ROS parameter server, connecting to the robot and setting up interfaces
         * for the ros_control framework.
         *
         * @tparam example of template param
         * 
         * @param[in] param example of function input param
         * 
         * @note example of note
         *
         * @return True if successful, false otherwise.
         */
        EigenEigenSolver();

        virtual ~EigenEigenSolver();

        /// Set solver tolerance
        void SetTol(double tol);

        /// Set maximum number of iterations
        void SetMaxIter(int max_iter);
        /// Set the number of required eigenmodes
        void SetNumModes(int num_eigs);
        /// Set operator for standard eigenvalue problem
        void SetOperator(const Eigen::MatrixXd& op);
        /// Set operator for generalized eigenvalue problem
        void SetOperators(const Eigen::MatrixXd& op, const Eigen::MatrixXd& opB);

        /// Customize object with options set
        void Customize(bool customize = true) const;

        /// Solve the eigenvalue problem for the specified number of eigenvalues
        void Solve();

        /// Get the number of converged eigenvalues
        int GetNumConverged();

        /// Get the corresponding eigenvalue
        void GetEigenvalue(unsigned int i, double& lr) const;
        void GetEigenvalue(unsigned int i, double& lr, double& lc) const;

        /// Get the corresponding eigenvector
        void GetEigenvector(unsigned int i, mfem::Vector& vr) const;
        void GetEigenvector(unsigned int i, mfem::Vector& vr, mfem::Vector& vc) const;

    protected:
    };
} // namespace gp_manifold

#endif // GPMANIFOLD_EIGEN_HPP