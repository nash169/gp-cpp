#include <Eigen/Core>
#include <iostream>
#include <sstream>

#include <magnum_dynamics/MagnumApp.hpp>
#include <utils_cpp/UtilsCpp.hpp>

int main(int argc, char** argv)
{
    magnum_dynamics::MagnumApp app({argc, argv});

    // Load mesh
    utils_cpp::FileManager io_manager;
    Eigen::MatrixXd vertices = io_manager.setFile("rsc/sol/sphere_mesh.mesh").read<Eigen::MatrixXd>("vertices", 3),
                    indices = io_manager.read<Eigen::MatrixXd>("elements", 2);

    size_t num_eig = 5;
    Eigen::MatrixXd eigenvectors(vertices.rows(), num_eig);

    for (size_t i = 0; i < num_eig; i++) {
        std::stringstream file_path;
        file_path << "rsc/sol/sphere_mode_0" << i;
        eigenvectors.col(i) = io_manager.setFile(file_path.str()).read<Eigen::MatrixXd>();
    }

    app.plot(vertices, eigenvectors.col(1), indices.block(0, 2, indices.rows(), 3))
        .setTransformation(Matrix4::translation({0.0f, 0.0f, 0.0f}));

    return app.exec();
}