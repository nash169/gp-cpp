#ifndef GP_MANIFOLD_SLEPC_SOLVER
#define GP_MANIFOLD_SLEPC_SOLVER

#include <iostream>
#include <memory>

#include <Eigen/Sparse>
#include <slepceps.h>

#include "gp_manifold/PetscMatrix.hpp"
#include "gp_manifold/PetscVector.hpp"

namespace gp_manifold {
    class SlepcSolver {
    public:
        SlepcSolver(int argc, char** argv)
        {
            // Init Slepc. This initialize MPI as well
            SlepcInitialize(&argc, &argv, (char*)0, NULL);

            // Create eigenvalue solver
            _ierr = EPSCreate(PETSC_COMM_WORLD, &_eps);
            CHKERRABORT(PETSC_COMM_WORLD, _ierr);
        }

        SlepcSolver& operators(std::unique_ptr<PetscMatrix> A, std::unique_ptr<PetscMatrix> B = nullptr)
        {
            _A = std::move(A);

            if (B)
                _B = std::move(B);
            else
                _B = std::make_unique<PetscMatrix>();

            _ierr = EPSSetOperators(_eps, *_A, *_B); //
            CHKERRABORT(PETSC_COMM_WORLD, _ierr);

            // Create Vector to store the solution
            _ierr = MatCreateVecs(*_A, NULL, &_xr);
            CHKERRABORT(PETSC_COMM_WORLD, _ierr);
            _ierr = MatCreateVecs(*_A, NULL, &_xi);
            CHKERRABORT(PETSC_COMM_WORLD, _ierr);

            return *this;
        }

        SlepcSolver& problem(const EPSProblemType& problem)
        {
            /* Set problem type */
            _ierr = EPSSetProblemType(_eps, problem);
            CHKERRABORT(PETSC_COMM_WORLD, _ierr);

            return *this;
        }

        SlepcSolver& spectrum(EPSWhich which)
        {
            _ierr = EPSSetWhichEigenpairs(_eps, which);
            CHKERRABORT(PETSC_COMM_WORLD, _ierr);

            return *this;
        }

        SlepcSolver& options(const std::string& prefix)
        {
            _ierr = EPSSetOptionsPrefix(_eps, prefix.c_str());
            CHKERRABORT(PETSC_COMM_WORLD, _ierr);

            return *this;
        }

        SlepcSolver& tolerance(const double& tol)
        {
            PetscInt max_its;

            _ierr = EPSGetTolerances(_eps, NULL, &max_its);
            CHKERRABORT(PETSC_COMM_WORLD, _ierr);
            // Work around uninitialized maximum iterations
            if (max_its == 0) {
                max_its = PETSC_DECIDE;
            }
            _ierr = EPSSetTolerances(_eps, tol, max_its);
            CHKERRABORT(PETSC_COMM_WORLD, _ierr);

            return *this;
        }

        SlepcSolver& iterations(const int& max_its)
        {
            double tol;

            _ierr = EPSGetTolerances(_eps, &tol, NULL);
            CHKERRABORT(PETSC_COMM_WORLD, _ierr);
            _ierr = EPSSetTolerances(_eps, tol, max_its);
            CHKERRABORT(PETSC_COMM_WORLD, _ierr);

            return *this;
        }

        SlepcSolver& modes(const int& num_eigs)
        {
            _ierr = EPSSetDimensions(_eps, num_eigs, PETSC_DECIDE, PETSC_DECIDE);
            CHKERRABORT(PETSC_COMM_WORLD, _ierr);

            return *this;
        }

        int converged()
        {
            PetscInt num_conv;
            _ierr = EPSGetConverged(_eps, &num_conv);
            // CHKERRABORT(PETSC_COMM_WORLD, _ierr);

            return static_cast<int>(num_conv);
        }

        void solve()
        {
            /* Set solver parameters at runtime */
            _ierr = EPSSetFromOptions(_eps);
            CHKERRABORT(PETSC_COMM_WORLD, _ierr);

            /* Solve the eigensystem */
            _ierr = EPSSolve(_eps);
            CHKERRABORT(PETSC_COMM_WORLD, _ierr);
        }

        PetscScalar eigenvalue(const size_t& i)
        {
            _ierr = EPSGetEigenvalue(_eps, i, &_kr, &_ki);
            // CHKERRQ(_ierr);

            return _kr;
        }

        Vec eigenvector(unsigned int i)
        {
            _ierr = EPSGetEigenvector(_eps, i, _xr, _xi);
            // CHKERRQ(_ierr);

            return _xr;
        }

        ~SlepcSolver()
        {
            _A.reset();

            if (_B)
                _B.reset();

            _ierr = EPSDestroy(&_eps);
            CHKERRABORT(PETSC_COMM_WORLD, _ierr);
            _ierr = VecDestroy(&_xr);
            CHKERRABORT(PETSC_COMM_WORLD, _ierr);
            _ierr = VecDestroy(&_xi);
            CHKERRABORT(PETSC_COMM_WORLD, _ierr);

            SlepcFinalize();
        }

    protected:
        /* Eigenproblem solver context */
        EPS _eps;
        EPSType _solver;
        EPSProblemType _problem;

        // Operators
        std::unique_ptr<PetscMatrix> _A, _B;

        // Error, tolerance, real and imaginary numbers
        PetscReal error, tol, re, im;

        // Real and imaginary eigenvalues
        PetscScalar _kr, _ki;

        // Real and imaginary eigenvectors
        Vec _xr, _xi;

        // Error code
        PetscErrorCode _ierr;
    };
} // namespace gp_manifold

#endif // GP_MANIFOLD_SLEPC_SOLVER