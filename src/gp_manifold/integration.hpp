#ifndef GPMANIFOLD_INTEGRATION_HPP
#define GPMANIFOLD_INTEGRATION_HPP

#include <Eigen/Sparse>
#include <mfem.hpp>
#include <vector>

/** @brief Eigen template specialization for vector conversion */
template <typename T>
struct VectorConverter {
    static mfem::Vector from(const Eigen::Matrix<T, Eigen::Dynamic, 1>& other)
    {
        mfem::Vector v(other.rows());

        for (size_t i = 0; i < v.Size(); i++)
            v(i) = other(i);

        return std::move(v);
    }

    static Eigen::Matrix<T, Eigen::Dynamic, 1> to(const mfem::Vector& other)
    {
        Eigen::Matrix<T, Eigen::Dynamic, 1> v(other.Size());

        for (size_t i = 0; i < v.Size(); i++)
            v(i) = other(i);

        return std::move(v);
    }
};

/** @brief Eigen template specialization for dense matrix conversion */
template <typename T>
struct DenseMatrixConverter {
    static mfem::DenseMatrix from(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& other)
    {
        mfem::DenseMatrix mat(other.rows(), other.cols());

        for (size_t j = 0; j < mat.Width(); j++)
            for (size_t i = 0; i < mat.Height(); i++)
                mat(i, j) = other(i, j);

        return mat;
    }

    static Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> to(const mfem::DenseMatrix& other)
    {
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat(other.Height(), other.Width());

        for (size_t j = 0; j < mat.cols(); j++)
            for (size_t i = 0; i < mat.rows(); i++)
                mat(i, j) = other(i, j);

        return mat;
    }
};

/** @brief Eigen template specialization for sparse matrix conversion */
template <class T>
struct SparseMatrixConverter {
    static mfem::SparseMatrix from(const Eigen::SparseMatrix<T, Eigen::RowMajor>& other)
    {
        return mfem::SparseMatrix(other.outerIndexPtr(), other.innerIndexPtr(), other.valuePtr(), other.rows(), other.cols());
    }

    static Eigen::SparseMatrix<T, Eigen::RowMajor> to(const mfem::SparseMatrix& other)
    {
        // MFEM memory info
        const int *I = other.GetI(), *J = other.GetJ();
        const T* Data = other.GetData();

        // Eigen triplet
        std::vector<Eigen::Triplet<double>> tripletList;
        tripletList.reserve(other.GetMemoryData().Capacity());

        for (size_t i = 0; i < other.Size(); i++) {
            for (size_t k = I[i], end = I[i + 1]; k < end; k++)
                tripletList.push_back(Eigen::Triplet<double>(i, J[k], Data[k]));
        }

        // Create Eigen sparse matrix
        Eigen::SparseMatrix<T, Eigen::RowMajor> mat(other.Height(), other.Width());
        mat.setFromTriplets(tripletList.begin(), tripletList.end());

        return mat;
    }
};

#endif // GP_MANIFOLD_INTEGRATION_HPP