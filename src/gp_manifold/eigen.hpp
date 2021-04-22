#ifndef GPMANIFOLD_EIGEN_HPP
#define GPMANIFOLD_EIGEN_HPP

#include <iostream>
#include <memory>
#include <vector>

#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>

#include <mfem.hpp>

#include "gp_manifold/integration.hpp"

namespace gp_manifold {
    class EigenEigenSolver {
    public:
        EigenEigenSolver() {}

        virtual ~EigenEigenSolver()
        {
            delete _A, _B, _S, _G;
        }

        /// Set solver tolerance
        void SetTol(double tol) {}

        /// Set maximum number of iterations
        void SetMaxIter(int max_iter) {}

        /// Set the number of required eigenmodes
        void SetNumModes(int num_eigs)
        {
        }

        /// Set operator for standard eigenvalue problem (A*x = lambda*x)
        void SetOperator(const mfem::Operator& A)
        {
            _A = new Eigen::SparseMatrix<double, Eigen::RowMajor>(SparseMatrixConverter<double>::to(static_cast<const mfem::BilinearForm&>(A).SpMat()));
        }

        /// Set operator for generalized eigenvalue problem (A*x = lambda*B*x)
        void SetOperators(const mfem::Operator& A, const mfem::Operator& B)
        {
            _A = new Eigen::SparseMatrix<double, Eigen::RowMajor>(SparseMatrixConverter<double>::to(static_cast<const mfem::BilinearForm&>(A).SpMat()));

            _B = new Eigen::SparseMatrix<double, Eigen::RowMajor>(SparseMatrixConverter<double>::to(static_cast<const mfem::BilinearForm&>(B).SpMat()));
        }

        /// Customize object with options set
        void Customize(bool customize = true) const {}

        /// Solve the eigenvalue problem for the specified number of eigenvalues
        void Solve()
        {
            if (!_B)
                _S = new Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(Eigen::MatrixXd(*_A));
            else
                _G = new Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd>(Eigen::MatrixXd(*_A), Eigen::MatrixXd(*_B));
        }

        /// Get the number of converged eigenvalues
        int GetNumConverged()
        {
            return 0;
        }

        /// Get the corresponding eigenvalue
        double GetEigenvalue(unsigned int i) const
        {
            if (!_B)
                return _S->eigenvalues()[i];
            else
                return _G->eigenvalues()[i];
        }

        Eigen::VectorXd GetEigenvalues() const
        {
            if (!_B)
                return _S->eigenvalues();
            else
                return _G->eigenvalues();
        }

        /// Get the corresponding eigenvector
        Eigen::VectorXd GetEigenvector(unsigned int i) const
        {
            if (!_B)
                return _S->eigenvectors().col(i);
            else
                return _G->eigenvectors().col(i);
        }

        Eigen::MatrixXd GetEigenvectors() const
        {
            if (!_B)
                return _S->eigenvectors();
            else
                return _G->eigenvectors();
        }

    protected:
        // Operators
        Eigen::SparseMatrix<double, Eigen::RowMajor>*_A = nullptr, *_B = nullptr;

        // Eigenvalue solution
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>* _S = nullptr;
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd>* _G = nullptr;
    };
} // namespace gp_manifold

#endif // GPMANIFOLD_EIGEN_HPP