#include <Eigen/Core>
#include <iostream>
#include <sstream>

#include <magnum_dynamics/MagnumApp.hpp>
#include <utils_cpp/UtilsCpp.hpp>

int main(int argc, char** argv)
{
    std::string mesh_name = "sphere";
    magnum_dynamics::MagnumApp app({argc, argv});
    utils_cpp::FileManager io_manager;

    // Load mesh
    Eigen::MatrixXd vertices = io_manager.setFile("rsc/modes/" + mesh_name + "_mesh.000000").read<Eigen::MatrixXd>("vertices", 3),
                    indices = io_manager.read<Eigen::MatrixXd>("elements", 2);

    // Load eigenvector
    size_t max_eig = 5;
    Eigen::VectorXd eigenvector(vertices.rows());

    std::stringstream file_path;
    file_path << "rsc/modes/" << mesh_name << "_mode_" << (argc >= max_eig ? 1 : argc) << ".000000";
    eigenvector = io_manager.setFile(file_path.str()).read<Eigen::MatrixXd>("", 5);

    // Plot eigenvector
    app.plot(vertices, eigenvector, indices.block(0, 2, indices.rows(), 3))
        .setTransformation(Matrix4::translation({0.0f, 0.0f, 0.0f}));

    return app.exec();
}