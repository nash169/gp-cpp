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
    // Eigen::MatrixXd vertices = io_manager.setFile("rsc/mesh/batman.msh").read<Eigen::MatrixXd>("$Nodes", 2, "$EndNodes"),
    //                 indices = io_manager.read<Eigen::MatrixXd>("$Elements", 2, "$EndElements").array() - 1;
    // Eigen::MatrixXd vertices = io_manager.setFile("rsc/sol/batman_mesh.mesh").read<Eigen::MatrixXd>("vertices", 3),
    //                 indices = io_manager.read<Eigen::MatrixXd>("elements", 2);

    Eigen::MatrixXd vertices = io_manager.setFile("rsc/data/mesh_vertices.csv").read<Eigen::MatrixXd>(),
                    indices = io_manager.setFile("rsc/data/mesh_faces.csv").read<Eigen::MatrixXd>();

    // size_t num_eig = 2;
    // Eigen::MatrixXd eigenvectors(vertices.rows(), num_eig);

    // for (size_t i = 0; i < num_eig; i++) {
    //     std::stringstream file_path;
    //     file_path << "rsc/sol/batman_mode_0" << i;
    //     eigenvectors.col(i) = io_manager.setFile(file_path.str()).read<Eigen::MatrixXd>();
    // }

    Eigen::VectorXd ground_truth = io_manager.setFile("rsc/data/ground_truth.csv").read<Eigen::MatrixXd>();

    // app.plot(vertices.block(0, 1, vertices.rows(), 3), eigenvectors.col(1), indices.block(0, 5, indices.rows(), 3))
    //     .setTransformation(Matrix4::translation({0.0f, 0.0f, 0.0f}));
    // app.plot(vertices, ground_truth, indices.block(0, 2, indices.rows(), 3)) // eigenvectors.col(1)  ground_truth
    //     .setTransformation(Matrix4::translation({0.0f, 0.0f, 0.0f}));

    app.plot(vertices, ground_truth, indices)
        .setTransformation(Matrix4::scaling({0.05, 0.05, 0.05}));

    // Eigen::Vector3f center = vertices.colwise().mean().cast<float>() * 0.05;
    // (*app.camera())
    //     .setCenter(Vector3(center))
    //     .setPose({10., 0., 5.});

    return app.exec();
}