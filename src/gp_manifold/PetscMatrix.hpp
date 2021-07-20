#ifndef GP_MANIFOLD_PETSC_MATRIX
#define GP_MANIFOLD_PETSC_MATRIX

#include <iostream>

#include <Eigen/Sparse>
#include <petsc.h>

// // forward declarations of PETSc internal structs
// struct _p_Vec;
// struct _p_Mat;

namespace gp_manifold {
    // namespace petsc {
    //     typedef struct ::_p_Vec* Vec;
    //     typedef struct ::_p_Mat* Mat;
    // } // namespace petsc

    class PetscMatrix {
    public:
        PetscMatrix()
        {
            A = nullptr;
        }

        // PetscMatrix()
        // {

        // }

        PetscMatrix(const Eigen::SparseMatrix<double, Eigen::RowMajor>& S)
        {
            ConvertOperator(PETSC_COMM_SELF, S, &A);
        }

        ~PetscMatrix()
        {
            PetscErrorCode ierr;
            ierr = MatDestroy(&A);
            CHKERRV(ierr);
        }

        Mat get()
        {
            return A;
        }

        operator Mat() const { return A; }

    protected:
        Mat A;

        int ConvertOperator(MPI_Comm comm, const Eigen::SparseMatrix<double, Eigen::RowMajor>& S, Mat* A)
        {
            /* from SparseMatrix to SEQAIJ -> always pass through host for now */
            Mat B;
            PetscScalar* pdata;
            PetscInt *pii, *pjj, *oii;
            PetscMPIInt size;
            PetscErrorCode ierr;

            int m = S.rows();
            int n = S.cols();
            const int* ii = S.outerIndexPtr();
            const int* jj = S.innerIndexPtr();
            const double* data = S.valuePtr();

            ierr = PetscMalloc1(m + 1, &pii);
            CHKERRQ(ierr);
            ierr = PetscMalloc1(ii[m], &pjj);
            CHKERRQ(ierr);
            ierr = PetscMalloc1(ii[m], &pdata);
            CHKERRQ(ierr);

            pii[0] = ii[0];
            for (int i = 0; i < m; i++) {
                bool issorted = true;
                pii[i + 1] = ii[i + 1];
                for (int j = ii[i]; j < ii[i + 1]; j++) {
                    pjj[j] = jj[j];
                    if (issorted && j != ii[i]) {
                        issorted = (pjj[j] > pjj[j - 1]);
                    }
                    pdata[j] = data[j];
                }
                if (!issorted) {
                    ierr = PetscSortIntWithScalarArray(pii[i + 1] - pii[i], pjj + pii[i], pdata + pii[i]);
                    CHKERRQ(ierr);
                }
            }
            ierr = MPI_Comm_size(comm, &size);
            CHKERRQ(ierr);

            if (size == 1) {
                ierr = MatCreateSeqAIJWithArrays(comm, m, n, pii, pjj, pdata, &B);
                CHKERRQ(ierr);
                oii = NULL;
            }
            else // block diagonal constructor
            {
                ierr = PetscCalloc1(m + 1, &oii);
                CHKERRQ(ierr);
                ierr = MatCreateMPIAIJWithSplitArrays(comm, m, n, PETSC_DECIDE,
                    PETSC_DECIDE,
                    pii, pjj, pdata, oii, NULL, NULL, &B);
                CHKERRQ(ierr);
            }
            // void* ptrs[4] = {pii, pjj, pdata, oii};
            // const char* names[4] = {"_mfem_csr_pii",
            //     "_mfem_csr_pjj",
            //     "_mfem_csr_pdata",
            //     "_mfem_csr_oii"};
            // for (int i = 0; i < 4; i++) {
            //     PetscContainer c;

            //     ierr = PetscContainerCreate(PETSC_COMM_SELF, &c);
            //     CHKERRQ(ierr);
            //     ierr = PetscContainerSetPointer(c, ptrs[i]);
            //     CHKERRQ(ierr);
            //     ierr = PetscContainerSetUserDestroy(c, __mfem_array_container_destroy);
            //     CHKERRQ(ierr);
            //     ierr = PetscObjectCompose((PetscObject)(B), names[i], (PetscObject)c);
            //     CHKERRQ(ierr);
            //     ierr = PetscContainerDestroy(&c);
            //     CHKERRQ(ierr);
            // }

            *A = B;

            return 1;
        }
    };
} // namespace gp_manifold

#endif // GP_MANIFOLD_PETSC_MATRIX