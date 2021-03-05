#include "eigen.hpp"

namespace gp_manifold {
    EigenEigenSolver::EigenEigenSolver() {}

    EigenEigenSolver::~EigenEigenSolver() {}

    void EigenEigenSolver::SetTol(double tol) {}

    void EigenEigenSolver::SetMaxIter(int max_iter) {}

    void EigenEigenSolver::SetNumModes(int num_eigs) {}

    void EigenEigenSolver::SetOperator(const Eigen::MatrixXd& op) {}

    void EigenEigenSolver::SetOperators(const Eigen::MatrixXd& op, const Eigen::MatrixXd& opB) {}

    void EigenEigenSolver::Customize(bool customize) const {}

    void EigenEigenSolver::Solve() {}

    int EigenEigenSolver::GetNumConverged()
    {
        return 1;
    }

    void EigenEigenSolver::GetEigenvalue(unsigned int i, double& lr) const {}

    void EigenEigenSolver::GetEigenvalue(unsigned int i, double& lr, double& lc) const {}

    void EigenEigenSolver::GetEigenvector(unsigned int i, mfem::Vector& vr) const {}

    void EigenEigenSolver::GetEigenvector(unsigned int i, mfem::Vector& vr, mfem::Vector& vc) const {}
} // namespace gp_manifold