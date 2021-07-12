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
    Eigen::MatrixXd nodes = io_manager.setFile("rsc/data/mesh_vertices.csv").read<Eigen::MatrixXd>(),
                    indices = io_manager.setFile("rsc/data/mesh_faces.csv").read<Eigen::MatrixXd>();

    Eigen::VectorXd ground_truth = io_manager.setFile("rsc/data/ground_truth.csv").read<Eigen::MatrixXd>(),
                    gp_sol = io_manager.setFile("rsc/data/gp_sol.csv").read<Eigen::MatrixXd>(),
                    rgp_sol = io_manager.setFile("rsc/data/rgp_sol.csv").read<Eigen::MatrixXd>();

    std::cout << ground_truth.minCoeff() << " : " << ground_truth.maxCoeff() << std::endl;
    std::cout << gp_sol.minCoeff() << " : " << gp_sol.maxCoeff() << std::endl;
    std::cout << rgp_sol.minCoeff() << " : " << rgp_sol.maxCoeff() << std::endl;

    app.plot(nodes, ground_truth, indices).setTransformation(Matrix4::translation({0.0f, 0.0f, 0.0f}));
    app.plot(nodes, gp_sol, indices).setTransformation(Matrix4::translation({0.0f, -3.0f, 0.0f}));
    app.plot(nodes, rgp_sol, indices).setTransformation(Matrix4::translation({0.0f, 3.0f, 0.0f}));

    return app.exec();
}