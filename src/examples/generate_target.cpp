#include <Eigen/Core>
#include <iostream>
#include <random>
#include <sstream>
#include <utils_lib/FileManager.hpp>

using namespace utils_lib;

int main(int argc, char** argv)
{
    constexpr int RAND_NUMS_TO_GENERATE = 100;

    std::string mesh_name = (argc > 1) ? argv[1] : "sphere";

    FileManager io_manager;

    Eigen::MatrixXd nodes = io_manager.setFile("rsc/truth/" + mesh_name + "_vertices.csv").read<Eigen::MatrixXd>(),
                    reference(RAND_NUMS_TO_GENERATE, nodes.cols());

    Eigen::VectorXd ground_truth = io_manager.setFile("rsc/truth/" + mesh_name + "_truth.csv").read<Eigen::MatrixXd>(),
                    target(RAND_NUMS_TO_GENERATE);

    // Extract reference and ground truth
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_int_distribution<int> distr(0, nodes.rows());

    for (size_t i = 0; i < RAND_NUMS_TO_GENERATE; i++) {
        int index = distr(eng);
        reference.row(i) = nodes.row(index);
        target(i) = ground_truth(index);
    }

    // Save target
    io_manager.setFile("rsc/truth/" + mesh_name + "_reference.csv").write(reference);
    io_manager.setFile("rsc/truth/" + mesh_name + "_target.csv").write(target);

    return 0;
}