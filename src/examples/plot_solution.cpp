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
                    gp = io_manager.setFile("rsc/solutions/" + mesh_name + "_gp.csv").read<Eigen::MatrixXd>(),
                    rgp_fem = io_manager.setFile("rsc/solutions/" + mesh_name + "_rgp_fem.csv").read<Eigen::MatrixXd>(),
                    rgp_diffusion = io_manager.setFile("rsc/solutions/" + mesh_name + "_rgp_diffusion.csv").read<Eigen::MatrixXd>();

    // std::cout << ground_truth.minCoeff() << " : " << ground_truth.maxCoeff() << std::endl;
    // std::cout << gp_sol.minCoeff() << " : " << gp_sol.maxCoeff() << std::endl;
    // std::cout << rgp_sol.minCoeff() << " : " << rgp_sol.maxCoeff() << std::endl;

    app.plot(nodes, ground_truth, faces).setTransformation(Matrix4::translation({0.0f, 0.0f, 1.5f}));
    app.plot(nodes, gp, faces).setTransformation(Matrix4::translation({0.0f, -3.0f, -1.5f}));
    app.plot(nodes, rgp_fem, faces).setTransformation(Matrix4::translation({0.0f, 0.0f, -1.5f}));
    app.plot(nodes, rgp_diffusion, faces).setTransformation(Matrix4::translation({0.0f, 3.0f, -1.5f}));

    return app.exec();
}