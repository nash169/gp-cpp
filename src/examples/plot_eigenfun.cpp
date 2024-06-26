#include <Eigen/Core>
#include <iostream>
#include <sstream>

#include <graphics_lib/Graphics.hpp>
#include <utils_lib/FileManager.hpp>

using namespace utils_lib;

int main(int argc, char** argv)
{
    std::string mesh_name = (argc > 1) ? argv[1] : "sphere";
    graphics_lib::Graphics app({argc, argv});
    FileManager io_manager;

    // Load mesh
    Eigen::MatrixXd fem_vertices = io_manager.setFile("outputs/modes/fem_" + mesh_name + "_mesh.000000").read<Eigen::MatrixXd>("vertices", 3),
                    fem_indices = io_manager.read<Eigen::MatrixXd>("elements", 2),
                    diffusion_vertices = io_manager.setFile("rsc/" + mesh_name + ".msh").read<Eigen::MatrixXd>("$Nodes", 2, "$EndNodes"),
                    diffusion_indices = io_manager.read<Eigen::MatrixXd>("$Elements", 2, "$EndElements").array() - 1;

    fem_indices.block(0, 0, fem_indices.rows(), 3) = fem_indices.block(0, 2, fem_indices.rows(), 3);
    fem_indices.conservativeResize(fem_indices.rows(), 3);

    diffusion_vertices.block(0, 0, diffusion_vertices.rows(), 3) = diffusion_vertices.block(0, 1, diffusion_vertices.rows(), 3);
    diffusion_vertices.conservativeResize(diffusion_vertices.rows(), 3);

    diffusion_indices.block(0, 0, diffusion_indices.rows(), 3) = diffusion_indices.block(0, 5, diffusion_indices.rows(), 3);
    diffusion_indices.conservativeResize(diffusion_indices.rows(), 3);

    // Load eigenvector
    Eigen::MatrixXd fem_vecs = io_manager.setFile("outputs/modes/fem_" + mesh_name + "_modes.000000").read<Eigen::MatrixXd>("modes", 2),
                    diffusion_vecs = io_manager.setFile("outputs/modes/diffusion_" + mesh_name + "_modes.000000").read<Eigen::MatrixXd>("modes", 2);

    // Plot eigenvector
    int fun_to_plot = (argc > 2) ? std::stoi(argv[2]) : 1;

    double min = (fem_vecs.col(fun_to_plot).minCoeff() < diffusion_vecs.col(fun_to_plot).minCoeff()) ? fem_vecs.col(fun_to_plot).minCoeff() : diffusion_vecs.col(fun_to_plot).minCoeff(),
           max = (fem_vecs.col(fun_to_plot).maxCoeff() > diffusion_vecs.col(fun_to_plot).maxCoeff()) ? fem_vecs.col(fun_to_plot).maxCoeff() : diffusion_vecs.col(fun_to_plot).maxCoeff();

    std::cout << min << " - " << max << std::endl;

    min = -0.05; // min - min * 0.01;
    max = 0.05; // max + max * 0.01;

    app.surf(fem_vertices, fem_vecs.col(fun_to_plot), fem_indices, min, max)
        .setTransformation(Matrix4::translation({0.0f, 5.0f, 0.0f}));

    app.surf(diffusion_vertices, diffusion_vecs.col(fun_to_plot), diffusion_indices, min, max)
        .setTransformation(Matrix4::translation({0.0f, -5.0f, 0.0f}));

    return app.exec();
}