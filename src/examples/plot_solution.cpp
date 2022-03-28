#include <Eigen/Core>
#include <iostream>
#include <sstream>

#include <graphics_lib/Graphics.hpp>
#include <utils_lib/FileManager.hpp>

using namespace utils_lib;

int main(int argc, char** argv)
{
    // Load data
    FileManager io_manager;

    std::string mesh_name = (argc > 1) ? argv[1] : "sphere";
    Eigen::MatrixXd nodes = io_manager.setFile("rsc/truth/" + mesh_name + "_vertices.csv").read<Eigen::MatrixXd>(),
                    faces = io_manager.setFile("rsc/truth/" + mesh_name + "_faces.csv").read<Eigen::MatrixXd>();

    std::string type = (argc > 2) ? argv[2] : "diffusion";
    Eigen::VectorXd sol(nodes.rows());
    if (type.compare("truth") == 0)
        sol = io_manager.setFile("rsc/truth/" + mesh_name + "_truth.csv").read<Eigen::MatrixXd>();
    else if (type.compare("ambient") == 0)
        sol = io_manager.setFile("rsc/solutions/ambient_" + mesh_name + "_gp.csv").read<Eigen::MatrixXd>();
    else if (type.compare("fem") == 0)
        sol = io_manager.setFile("rsc/solutions/fem_" + mesh_name + "_rgp.csv").read<Eigen::MatrixXd>("sol", 2);
    else if (type.compare("diffusion") == 0)
        sol = io_manager.setFile("rsc/solutions/diffusion_" + mesh_name + "_rgp.csv").read<Eigen::MatrixXd>();

    // Set min and max value for coloring
    double min = sol.minCoeff() - sol.minCoeff() * 0.05,
           max = sol.maxCoeff() + sol.maxCoeff() * 0.05;
    std::cout << min << " - " << max << std::endl;

    // Graphics
    graphics_lib::Graphics app({argc, argv});
    app.setBackground("white").surf(nodes, sol, faces, min, max);

    return app.exec();
}