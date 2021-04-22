#include <Eigen/Sparse>
#include <iostream>
#include <mfem.hpp>

#include <gp_manifold/integration.hpp>

int main(int argc, char const* argv[])
{
    int I[] = {0, 2, 4, 6, 7}, J[] = {0, 1, 1, 3, 2, 4, 5}, m = 4, n = 6;
    double Data[] = {10, 20, 30, 40, 50, 70, 80};

    // mfem::Memory<int> I, J;
    // mfem::Memory<double> A;

    // I.Wrap(i, m + 1, true);
    // J.Wrap(j, I[m], true);
    // A.Wrap(data, I[m], true);

    // mfem::SparseMatrix* mat = new mfem::SparseMatrix(I, J, Data, m, n, false, false, true);
    mfem::SparseMatrix mat(I, J, Data, m, n, false, false, true);

    // std::cout << mat->Elem(2, 2) << std::endl;

    // std::vector<Eigen::Triplet<double>> tripletList;
    // tripletList.reserve(mat->GetMemoryData().Capacity());

    // std::cout << mat->GetMemoryData().Capacity() << std::endl;

    // for (size_t i = 0; i < mat->Size(); i++) {
    //     for (size_t k = I[i], end = I[i + 1]; k < end; k++) {
    //         std::cout << i << " " << J[k] << " " << Data[k] << std::endl;
    //         tripletList.push_back(Eigen::Triplet<double>(i, J[k], Data[k]));
    //     }
    // }

    // std::cout << tripletList.size() << std::endl;

    // std::cout << mat->Height() << " " << mat->Width() << std::endl;

    // Create Eigen sparse matrix
    // Eigen::SparseMatrix<double> mat_eigen(100, 100);

    // mat_eigen.setFromTriplets(tripletList.begin(), tripletList.end());

    // std::cout << "HELLO" << std::endl;

    Eigen::SparseMatrix<double, Eigen::RowMajor> mat_eigen = SparseMatrixConverter<double>::to(mat);
    // Eigen::SparseMatrix<double, Eigen::RowMajor> mat_eigen = MyFun<double>(*mat);
    // Eigen::SparseMatrix<double, Eigen::RowMajor> mat_eigen(Eigen::SparseMatrix<double, Eigen::RowMajor>(*mat));

    // mat_eigen.makeCompressed();

    // std::cout << mat_eigen.innerSize() << std::endl;

    // for (int k; k < mat_eigen.innerSize() + 1; k++) {
    //     // int j = mat_eigen.innerIndexPtr()[k];
    //     // double v = mat_eigen.valuePtr()[k];
    //     std::cout << mat_eigen.valuePtr()[k] << std::endl;
    //     // v is value of the element at position (j,i)
    // }

    // std::cout << mat_eigen.innerSize() << std::endl;

    Eigen::MatrixXd dense(mat_eigen);

    std::cout << dense << std::endl;

    // delete mat;

    return 0;
}
