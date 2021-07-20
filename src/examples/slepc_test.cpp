/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2021, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Standard symmetric eigenproblem corresponding to the Laplacian operator in 1 dimension.\n\n"
                     "The command line options are:\n"
                     "  -n <n>, where <n> = number of grid subdivisions = matrix dimension.\n\n";

#include <gp_manifold/PetscMatrix.hpp>
#include <iostream>
#include <kernel_lib/Kernel.hpp>
#include <slepceps.h>

using namespace gp_manifold;
using namespace kernel_lib;
// Kernel parameters
struct Params {
    struct kernel : public defaults::kernel {
        PARAM_SCALAR(double, sf, 0.405465);
        PARAM_SCALAR(double, sn, 0.693147);
    };

    struct exp_sq : public defaults::exp_sq {
        PARAM_SCALAR(double, l, -0.356675);
    };
};

int main(int argc, char** argv)
{
    int num_procs, myid;
    Mat A; /* problem matrix */
    EPS eps; /* eigenproblem solver context */
    EPSType type;
    PetscReal error, tol, re, im;
    PetscScalar kr, ki;
    Vec xr, xi;
    PetscInt n = 10, i, Istart, Iend, nev, maxit, its, nconv;
    PetscErrorCode ierr;

    // Init Slepc. This initialize MPI as well
    ierr = SlepcInitialize(&argc, &argv, (char*)0, NULL);
    if (ierr)
        return ierr;

    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    Eigen::SparseMatrix<double, Eigen::RowMajor> mat;

    // if (0 == myid) {
    constexpr int dim = 2, num_samples = 10;

    Eigen::MatrixXd X = Eigen::MatrixXd::Random(num_samples, dim);

    utils::Graph graph;

    mat = graph.kNearestWeighted(X, kernels::SquaredExp<Params>(), 3, 1);
    //

    PetscMatrix PetscMat(mat);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    PetscPrintf(PETSC_COMM_WORLD, "\nThe Kershaw matrix:\n\n");
    MatView(PetscMat.get(), PETSC_VIEWER_STDOUT_WORLD);

    PetscMat.~PetscMatrix();

    ierr = MatCreate(PETSC_COMM_WORLD, &A);
    CHKERRQ(ierr);
    ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n);
    CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);
    CHKERRQ(ierr);
    ierr = MatSetUp(A);
    CHKERRQ(ierr);

    ierr = MatGetOwnershipRange(A, &Istart, &Iend);
    CHKERRQ(ierr);
    for (i = Istart; i < Iend; i++) {
        if (i > 0) {
            ierr = MatSetValue(A, i, i - 1, -1.0, INSERT_VALUES);
            CHKERRQ(ierr);
        }
        if (i < n - 1) {
            ierr = MatSetValue(A, i, i + 1, -1.0, INSERT_VALUES);
            CHKERRQ(ierr);
        }
        ierr = MatSetValue(A, i, i, 2.0, INSERT_VALUES);
        CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    ierr = MatCreateVecs(A, NULL, &xr);
    CHKERRQ(ierr);
    ierr = MatCreateVecs(A, NULL, &xi);
    CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
     Create eigensolver context
  */
    ierr = EPSCreate(PETSC_COMM_WORLD, &eps);
    CHKERRQ(ierr);

    /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
    ierr = EPSSetOperators(eps, A, NULL);
    CHKERRQ(ierr);
    ierr = EPSSetProblemType(eps, EPS_HEP);
    CHKERRQ(ierr);

    /*
     Set solver parameters at runtime
  */
    ierr = EPSSetFromOptions(eps);
    CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    ierr = EPSSolve(eps);
    CHKERRQ(ierr);
    /*
     Optional: Get some information from the solver and display it
  */
    ierr = EPSGetIterationNumber(eps, &its);
    CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, " Number of iterations of the method: %D\n", its);
    CHKERRQ(ierr);
    ierr = EPSGetType(eps, &type);
    CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, " Solution method: %s\n\n", type);
    CHKERRQ(ierr);
    ierr = EPSGetDimensions(eps, &nev, NULL, NULL);
    CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, " Number of requested eigenvalues: %D\n", nev);
    CHKERRQ(ierr);
    ierr = EPSGetTolerances(eps, &tol, &maxit);
    CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, " Stopping condition: tol=%.4g, maxit=%D\n", (double)tol, maxit);
    CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
     Get number of converged approximate eigenpairs
  */
    ierr = EPSGetConverged(eps, &nconv);
    CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, " Number of converged eigenpairs: %D\n\n", nconv);
    CHKERRQ(ierr);

    if (nconv > 0) {
        /*
       Display eigenvalues and relative errors
    */
        ierr = PetscPrintf(PETSC_COMM_WORLD,
            "           k          ||Ax-kx||/||kx||\n"
            "   ----------------- ------------------\n");
        CHKERRQ(ierr);

        for (i = 0; i < nconv; i++) {
            /*
        Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
        ki (imaginary part)
      */
            ierr = EPSGetEigenpair(eps, i, &kr, &ki, xr, xi);
            CHKERRQ(ierr);
            /*
         Compute the relative error associated to each eigenpair
      */
            ierr = EPSComputeError(eps, i, EPS_ERROR_RELATIVE, &error);
            CHKERRQ(ierr);

#if defined(PETSC_USE_COMPLEX)
            re = PetscRealPart(kr);
            im = PetscImaginaryPart(kr);
#else
            re = kr;
            im = ki;
#endif
            if (im != 0.0) {
                ierr = PetscPrintf(PETSC_COMM_WORLD, " %9f%+9fi %12g\n", (double)re, (double)im, (double)error);
                CHKERRQ(ierr);
            }
            else {
                ierr = PetscPrintf(PETSC_COMM_WORLD, "   %12f       %12g\n", (double)re, (double)error);
                CHKERRQ(ierr);
            }
        }
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");
        CHKERRQ(ierr);
    }

    PetscPrintf(PETSC_COMM_WORLD, "\nEigenvector:\n\n");
    VecView(xr, PETSC_VIEWER_STDOUT_WORLD);

    /*
     Free work space
  */
    ierr = EPSDestroy(&eps);
    CHKERRQ(ierr);
    ierr = MatDestroy(&A);
    CHKERRQ(ierr);
    ierr = VecDestroy(&xr);
    CHKERRQ(ierr);
    ierr = VecDestroy(&xi);
    CHKERRQ(ierr);
    ierr = SlepcFinalize();
    return ierr;
}
