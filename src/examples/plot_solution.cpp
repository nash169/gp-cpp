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

    // Load ground truth, target and relative nodes
    Eigen::MatrixXd nodes = io_manager.setFile("rsc/truth/" + mesh_name + "_vertices.csv").read<Eigen::MatrixXd>(),
                    faces = io_manager.setFile("rsc/truth/" + mesh_name + "_faces.csv").read<Eigen::MatrixXd>();

    Eigen::VectorXd ground_truth = io_manager.setFile("rsc/truth/" + mesh_name + "_truth.csv").read<Eigen::MatrixXd>(),
                    gp = io_manager.setFile("rsc/solutions/ambient_" + mesh_name + "_gp.csv").read<Eigen::MatrixXd>(),
                    rgp_fem = io_manager.setFile("rsc/solutions/fem_" + mesh_name + "_rgp.csv").read<Eigen::MatrixXd>(),
                    rgp_diffusion = io_manager.setFile("rsc/solutions/diffusion_" + mesh_name + "_rgp.csv").read<Eigen::MatrixXd>();

    // Plot
    Eigen::Vector4d min_coeff(ground_truth.minCoeff(), gp.minCoeff(), rgp_fem.minCoeff(), rgp_diffusion.minCoeff()),
        max_coeff(ground_truth.maxCoeff(), gp.maxCoeff(), rgp_fem.maxCoeff(), rgp_diffusion.maxCoeff());

    double min = min_coeff.minCoeff() - min_coeff.minCoeff() * 0.05,
           max = max_coeff.maxCoeff() + max_coeff.maxCoeff() * 0.05;

    std::cout << min << " - " << max << std::endl;
    std::cout << ground_truth.minCoeff() << " - " << ground_truth.maxCoeff() << std::endl;
    std::cout << gp.minCoeff() << " - " << gp.maxCoeff() << std::endl;
    std::cout << rgp_fem.minCoeff() << " - " << rgp_fem.maxCoeff() << std::endl;
    std::cout << rgp_diffusion.minCoeff() << " - " << rgp_diffusion.maxCoeff() << std::endl;

    app.surf(nodes, ground_truth, faces, min, max)
        .setTransformation(Matrix4::translation(Vector3(0, 0, 2)));

    app.surf(nodes, gp, faces, min, max)
        .setTransformation(Matrix4::translation(Vector3(0, 0, -2)));

    app.surf(nodes, rgp_fem, faces, min, max)
        .setTransformation(Matrix4::translation(Vector3(0, -10, -2)));

    app.surf(nodes, rgp_diffusion, faces, min, max)
        .setTransformation(Matrix4::translation(Vector3(0, 10, -2)));

    return app.exec();
}