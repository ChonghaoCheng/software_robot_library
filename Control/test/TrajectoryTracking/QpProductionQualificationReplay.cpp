#include <Math/BoxAwareActiveSet.h>
#include <Math/QpAcceptance.h>

#include <Eigen/Dense>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

Eigen::MatrixXd read_csv(const std::filesystem::path &path)
{
    std::ifstream stream(path);
    if(!stream) throw std::runtime_error("Unable to open " + path.string());
    std::vector<std::vector<double>> rows;
    std::string line;
    while(std::getline(stream, line))
    {
        if(line.empty()) continue;
        std::vector<double> row;
        std::stringstream parser(line);
        std::string field;
        while(std::getline(parser, field, ',')) row.push_back(std::stod(field));
        if(!rows.empty() && row.size() != rows.front().size())
            throw std::runtime_error("Ragged CSV: " + path.string());
        rows.push_back(std::move(row));
    }
    if(rows.empty()) return Eigen::MatrixXd(0, 0);
    Eigen::MatrixXd result(rows.size(), rows.front().size());
    for(Eigen::Index row = 0; row < result.rows(); ++row)
        for(Eigen::Index column = 0; column < result.cols(); ++column)
            result(row, column) = rows[static_cast<std::size_t>(row)]
                                      [static_cast<std::size_t>(column)];
    return result;
}

void write_vector(const std::filesystem::path &path, const Eigen::VectorXd &value)
{
    std::ofstream stream(path);
    if(!stream) throw std::runtime_error("Unable to write " + path.string());
    stream << std::setprecision(17);
    for(Eigen::Index index = 0; index < value.size(); ++index)
        stream << value(index) << '\n';
}

double primal_violation(const Eigen::MatrixXd &B, const Eigen::VectorXd &z,
                        const Eigen::VectorXd &solution)
{
    if(B.rows() == 0) return 0.0;
    return std::max(0.0, (B * solution - z).maxCoeff());
}

} // namespace

int main(int argc, char **argv)
{
    if(argc != 2)
    {
        std::cerr << "Usage: qp_production_qualification_replay QP_DIR\n";
        return 2;
    }
    try
    {
        const std::filesystem::path directory(argv[1]);
        const Eigen::MatrixXd H = read_csv(directory / "H.csv");
        const Eigen::VectorXd f = read_csv(directory / "f.csv").col(0);
        const Eigen::MatrixXd B = read_csv(directory / "B.csv");
        const Eigen::VectorXd z = read_csv(directory / "z.csv").col(0);
        const Eigen::VectorXd seed = read_csv(directory / "seed.csv").col(0);

        SolverOptions<double> options;
        options.method = "active set";
        options.maxSteps = 50;
        options.stepSizeTolerance = 1e-5;
        const auto result = solve_box_aware_active_set(H, f, B, z, seed, options);
        const double violation = primal_violation(B, z, result.solution);
        RobotLibrary::Math::require_qp_result_accepted(
            result.solver, result.solution, violation, 1e-8,
            "QP production qualification replay");
        write_vector(directory / "qualified_box_aware_solution.csv", result.solution);

        std::cout << std::setprecision(17)
                  << "{\"termination_reason\":\"Converged\","
                  << "\"number_of_steps\":" << result.solver.numberOfSteps << ','
                  << "\"final_step_size\":" << result.solver.finalStepSize << ','
                  << "\"objective\":" << result.solver.objectiveFunction << ','
                  << "\"primal_violation\":" << violation << ','
                  << "\"active_set_changes\":" << result.solver.activeSetChanges << ','
                  << "\"maximum_active_set_size\":"
                  << result.solver.maximumActiveSetSize << ','
                  << "\"unique_active_constraints\":"
                  << result.solver.uniqueActiveConstraints << ','
                  << "\"cache_hits\":" << result.cacheHits << ','
                  << "\"cache_misses\":" << result.cacheMisses << ','
                  << "\"repeated_activations\":" << result.repeatedActivations
                  << "}\n";
        return 0;
    }
    catch(const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
