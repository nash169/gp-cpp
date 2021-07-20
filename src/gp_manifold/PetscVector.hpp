#ifndef GP_MANIFOLD_PETSC_VECTOR
#define GP_MANIFOLD_PETSC_VECTOR

#include <iostream>

#include <Eigen/Core>
#include <petsc.h>

namespace gp_manifold {
    class PetscVector {
    public:
        PetscVector()
        {
            _V = nullptr;
        }

        PetscVector(const Vec& V)
        {
            _V = V; // this has to be moved
        }

        PetscVector(const Eigen::Matrix<double, Eigen::Dynamic, 1>& V)
        {
            // make the conversion here
        }

        ~PetscVector()
        {
            _ierr = VecDestroy(&_V);
            CHKERRABORT(PETSC_COMM_WORLD, _ierr);
        }

        operator Vec() const { return _V; }

    protected:
        PetscErrorCode _ierr;

        Vec _V;
    };
} // namespace gp_manifold

#endif // GP_MANIFOLD_PETSC_VECTOR