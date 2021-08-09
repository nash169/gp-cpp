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
    Eigen::MatrixXd fem_vertices = io_manager.setFile("rsc/modes/fem_" + mesh_name + "_mesh.000000").read<Eigen::MatrixXd>("vertices", 3),
                    fem_indices = io_manager.read<Eigen::MatrixXd>("elements", 2),
                    diffusion_vertices = io_manager.setFile("rsc/meshes/" + mesh_name + ".msh").read<Eigen::MatrixXd>("$Nodes", 2, "$EndNodes"),
                    diffusion_indices = io_manager.read<Eigen::MatrixXd>("$Elements", 2, "$EndElements").array() - 1;

    fem_indices.block(0, 0, fem_indices.rows(), 3) = fem_indices.block(0, 2, fem_indices.rows(), 3);
    fem_indices.conservativeResize(fem_indices.rows(), 3);

    diffusion_vertices.block(0, 0, diffusion_vertices.rows(), 3) = diffusion_vertices.block(0, 1, diffusion_vertices.rows(), 3);
    diffusion_vertices.conservativeResize(diffusion_vertices.rows(), 3);

    diffusion_indices.block(0, 0, diffusion_indices.rows(), 3) = diffusion_indices.block(0, 5, diffusion_indices.rows(), 3);
    diffusion_indices.conservativeResize(diffusion_indices.rows(), 3);

    // Load eigenvector
    Eigen::MatrixXd fem_vecs = io_manager.setFile("rsc/modes/fem_" + mesh_name + "_modes.000000").read<Eigen::MatrixXd>("modes", 2),
                    diffusion_vecs = io_manager.setFile("rsc/modes/diffusion_" + mesh_name + "_modes.000000").read<Eigen::MatrixXd>("modes", 2);

    // Plot eigenvector
    int fun_to_plot = 4;

    // std::cout << fem_vecs.col(fun_to_plot).maxCoeff() << " - " << fem_vecs.col(fun_to_plot).minCoeff() << std::endl;
    // app.plot(fem_vertices, fem_vecs.col(fun_to_plot), fem_indices, -0.06, 0.06);

    std::cout << diffusion_vecs.col(fun_to_plot).maxCoeff() << " - " << diffusion_vecs.col(fun_to_plot).minCoeff() << std::endl;
    app.plot(diffusion_vertices, diffusion_vecs.col(fun_to_plot), diffusion_indices, -0.055, 0.055);

    return app.exec();
}