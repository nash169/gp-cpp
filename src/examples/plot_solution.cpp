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

    // Load ground truth, target and relative nodes
    Eigen::MatrixXd nodes = io_manager.setFile("rsc/truth/" + mesh_name + "_vertices.csv").read<Eigen::MatrixXd>(),
                    faces = io_manager.setFile("rsc/truth/" + mesh_name + "_faces.csv").read<Eigen::MatrixXd>();

    Eigen::VectorXd ground_truth = io_manager.setFile("rsc/truth/" + mesh_name + "_truth.csv").read<Eigen::MatrixXd>(),
                    gp = io_manager.setFile("rsc/solutions/ambient_" + mesh_name + "_gp.csv").read<Eigen::MatrixXd>(),
                    rgp_fem = io_manager.setFile("rsc/solutions/fem_" + mesh_name + "_rgp.csv").read<Eigen::MatrixXd>(),
                    rgp_diffusion = io_manager.setFile("rsc/solutions/diffusion_" + mesh_name + "_rgp.csv").read<Eigen::MatrixXd>();

    // Plot

    // std::cout << ground_truth.minCoeff() << " - " << ground_truth.maxCoeff() << std::endl;
    // app.plot(nodes, ground_truth, faces, -1.2, 1.2);

    // std::cout << gp.minCoeff() << " - " << gp.maxCoeff() << std::endl;
    // app.plot(nodes, gp, faces, -1.2, 1.2);

    // std::cout << rgp_fem.minCoeff() << " - " << rgp_fem.maxCoeff() << std::endl;
    // app.plot(nodes, rgp_fem, faces, -1.1, 1.1);

    std::cout << rgp_diffusion.minCoeff() << " - " << rgp_diffusion.maxCoeff() << std::endl;
    app.plot(nodes, rgp_diffusion, faces, -1.2, 1.2);

    return app.exec();
}