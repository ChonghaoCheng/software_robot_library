#include <Math/QPSolver.h>

#include <Eigen/Core>

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
    if(not stream)
    {
        throw std::runtime_error("Unable to open " + path.string());
    }
    std::vector<std::vector<double>> rows;
    std::string line;
    while(std::getline(stream, line))
    {
        if(line.empty())
        {
            continue;
        }
        std::vector<double> row;
        std::stringstream lineStream(line);
        std::string value;
        while(std::getline(lineStream, value, ','))
        {
            row.push_back(std::stod(value));
        }
        if(not rows.empty() and row.size() != rows.front().size())
        {
            throw std::runtime_error("Ragged CSV file " + path.string());
        }
        rows.push_back(std::move(row));
    }
    if(rows.empty())
    {
        return Eigen::MatrixXd(0, 0);
    }
    Eigen::MatrixXd result(
        static_cast<Eigen::Index>(rows.size()),
        static_cast<Eigen::Index>(rows.front().size()));
    for(Eigen::Index row = 0; row < result.rows(); ++row)
    {
        for(Eigen::Index column = 0; column < result.cols(); ++column)
        {
            result(row, column) = rows[static_cast<std::size_t>(row)]
                                      [static_cast<std::size_t>(column)];
        }
    }
    return result;
}

void append_row(Eigen::MatrixXd &matrix,
                Eigen::VectorXd &bounds,
                const Eigen::RowVectorXd &row,
                const double bound)
{
    const Eigen::Index oldRows = matrix.rows();
    matrix.conservativeResize(oldRows + 1, Eigen::NoChange);
    bounds.conservativeResize(oldRows + 1);
    matrix.row(oldRows) = row;
    bounds(oldRows) = bound;
}

void write_number(std::ostream &stream, const double value)
{
    if(std::isfinite(value))
    {
        stream << value;
    }
    else
    {
        stream << "null";
    }
}

} // namespace

int main(int argc, char **argv)
{
    if(argc != 3)
    {
        std::cerr << "Usage: n125_qp_active_set_replay SNAPSHOT_DIR V0|V1|V2\n";
        return 2;
    }
    try
    {
        const std::filesystem::path directory(argv[1]);
        const std::string variant(argv[2]);
        const Eigen::MatrixXd H = read_csv(directory / "H.csv");
        const Eigen::VectorXd f = read_csv(directory / "f.csv").col(0);
        const Eigen::MatrixXd originalB = read_csv(directory / "Bineq.csv");
        const Eigen::VectorXd originalZ = read_csv(directory / "zineq.csv").col(0);
        const Eigen::VectorXd nominal = read_csv(directory / "zNominal.csv").col(0);
        const Eigen::VectorXd fixedRates =
            read_csv(directory / "fixed_progress_rates.csv").col(0);
        const Eigen::Index dim = H.rows();
        const Eigen::Index horizon = fixedRates.size();
        const Eigen::Index progressOffset = dim - horizon;

        Eigen::MatrixXd Aeq(0, dim);
        Eigen::VectorXd yeq(0);
        Eigen::MatrixXd Bineq(0, dim);
        Eigen::VectorXd zineq(0);
        if(variant == "V0")
        {
            Bineq = originalB;
            zineq = originalZ;
        }
        else if(variant == "V1" or variant == "V2")
        {
            Aeq = Eigen::MatrixXd::Zero(horizon, dim);
            yeq = fixedRates;
            for(Eigen::Index stage = 0; stage < horizon; ++stage)
            {
                Aeq(stage, progressOffset + stage) = 1.0;
            }
            for(Eigen::Index row = 0; row < progressOffset; ++row)
            {
                append_row(Bineq, zineq, originalB.row(row), originalZ(row));
            }
            for(Eigen::Index row = dim; row < dim + progressOffset; ++row)
            {
                append_row(Bineq, zineq, originalB.row(row), originalZ(row));
            }
            if(variant == "V1")
            {
                append_row(Bineq, zineq,
                           originalB.row(2 * dim), originalZ(2 * dim));
                append_row(Bineq, zineq,
                           originalB.row(2 * dim + 1), originalZ(2 * dim + 1));
            }
        }
        else
        {
            throw std::invalid_argument("Unknown representation " + variant);
        }

        SolverOptions<double> options;
        options.method = "active set";
        options.stepSizeTolerance = 1e-5;
        options.maxSteps = 50;
        QPSolver<double> solver(options);
        const Eigen::VectorXd solution =
            solver.solve(H, f, Aeq, yeq, Bineq, zineq, nominal);
        const auto results = solver.results();
        const double inequalityViolation = Bineq.rows() == 0
            ? 0.0 : (Bineq * solution - zineq).maxCoeff();
        const double equalityResidual = Aeq.rows() == 0
            ? 0.0 : (Aeq * solution - yeq).cwiseAbs().maxCoeff();
        const double objective =
            0.5 * solution.dot(H * solution) + f.dot(solution);
        std::cout << std::setprecision(17)
                  << "{\"variant\":\"" << variant << "\","
                  << "\"status\":\"returned\","
                  << "\"number_of_steps\":" << results.numberOfSteps << ','
                  << "\"final_step_size\":";
        write_number(std::cout, results.finalStepSize);
        std::cout << ",\"max_inequality_violation\":";
        write_number(std::cout, inequalityViolation);
        std::cout << ",\"max_equality_residual\":";
        write_number(std::cout, equalityResidual);
        std::cout << ",\"objective\":";
        write_number(std::cout, objective);
        std::cout << ",\"solution_norm\":";
        write_number(std::cout, solution.norm());
        std::cout << "}\n";
        return 0;
    }
    catch(const std::exception &error)
    {
        std::cout << "{\"variant\":\"" << argv[2]
                  << "\",\"status\":\"exception\",\"message\":\""
                  << error.what() << "\"}\n";
        return 1;
    }
}
