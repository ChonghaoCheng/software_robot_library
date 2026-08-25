#include <Math/QPSolver.h>

#include <Eigen/Core>

#include <algorithm>
#include <chrono>
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
    if(not stream) throw std::runtime_error("Unable to open " + path.string());
    std::vector<std::vector<double>> rows;
    std::string line;
    while(std::getline(stream, line))
    {
        if(line.empty()) continue;
        std::vector<double> row;
        std::stringstream values(line);
        std::string value;
        while(std::getline(values, value, ',')) row.push_back(std::stod(value));
        if(not rows.empty() and row.size() != rows.front().size())
            throw std::runtime_error("Ragged CSV " + path.string());
        rows.push_back(std::move(row));
    }
    if(rows.empty()) return Eigen::MatrixXd(0, 0);
    Eigen::MatrixXd result(rows.size(), rows.front().size());
    for(Eigen::Index row = 0; row < result.rows(); ++row)
        for(Eigen::Index column = 0; column < result.cols(); ++column)
            result(row, column) = rows[static_cast<size_t>(row)]
                                      [static_cast<size_t>(column)];
    return result;
}

Eigen::VectorXd read_vector(const std::filesystem::path &path)
{
    const Eigen::MatrixXd matrix = read_csv(path);
    if(matrix.cols() != 1) throw std::runtime_error("Expected vector " + path.string());
    return matrix.col(0);
}

void write_vector(const std::filesystem::path &path, const Eigen::VectorXd &x)
{
    std::ofstream stream(path);
    stream << std::setprecision(17);
    for(Eigen::Index index = 0; index < x.size(); ++index) stream << x(index) << '\n';
}

std::vector<int> group_rows(const std::string &group, int dimension, int horizon)
{
    if(group == "ALL_ROWS")
    {
        std::vector<int> all;
        for(int row = 0; row < dimension; ++row) all.push_back(row);
        return all;
    }
    const int baseRows = 2 * (7 * horizon) + 1;
    std::vector<int> rows;
    for(int row = 0; row < baseRows; ++row) rows.push_back(row);
    if(group == "Q0") return rows;
    const int actionStart = baseRows + 6 * horizon;
    for(int stage = 0; stage < horizon; ++stage)
        for(int local = 0; local < (group == "Q1" ? 2 : 4); ++local)
            rows.push_back(actionStart + 4 * stage + local);
    if(group == "Q1" or group == "Q2") return rows;
    for(int stage = 0; stage < horizon; ++stage)
        for(int local = 0; local < 3; ++local)
            rows.push_back(baseRows + 6 * stage + local);
    if(group == "Q3") return rows;
    for(int stage = 0; stage < horizon; ++stage)
        for(int local = 3; local < 6; ++local)
            rows.push_back(baseRows + 6 * stage + local);
    if(group != "Q4" and group != "FULL"
       and group != "FULL_NO_NORMAL_GUARD"
       and group != "FULL_NO_SLEW"
       and group != "FULL_NO_MAGNITUDE")
        throw std::invalid_argument("Unknown group " + group);
    if(group == "FULL_NO_NORMAL_GUARD" or group == "FULL_NO_SLEW"
       or group == "FULL_NO_MAGNITUDE")
    {
        const int magnitudeBegin = actionStart;
        std::vector<int> filtered;
        filtered.reserve(rows.size());
        for(const int row : rows)
        {
            if(row < magnitudeBegin) { filtered.push_back(row); continue; }
            const int local = (row - magnitudeBegin) % 4;
            const bool magnitude = local < 2;
            const bool slew = local >= 2;
            if(group == "FULL_NO_NORMAL_GUARD") continue;
            if(group == "FULL_NO_SLEW" and slew) continue;
            if(group == "FULL_NO_MAGNITUDE" and magnitude) continue;
            filtered.push_back(row);
        }
        rows = std::move(filtered);
    }
    for(const int row : rows)
        if(row < 0 or row >= dimension)
            throw std::runtime_error("Constraint group row out of range");
    return rows;
}

Eigen::VectorXd projected_seed(Eigen::VectorXd seed,
                               const Eigen::MatrixXd &B,
                               const Eigen::VectorXd &z)
{
    for(int sweep = 0; sweep < 200; ++sweep)
    {
        double maximum = 0.0;
        for(Eigen::Index row = 0; row < B.rows(); ++row)
        {
            const double residual = B.row(row).dot(seed) - z(row);
            maximum = std::max(maximum, residual);
            const double normSquared = B.row(row).squaredNorm();
            if(residual > 1e-12 and normSquared > 1e-20)
                seed -= B.row(row).transpose() * ((residual + 1e-12) / normSquared);
        }
        if(maximum <= 1e-10) break;
    }
    return seed;
}

void polish(Eigen::VectorXd &x, const Eigen::MatrixXd &B,
            const Eigen::VectorXd &z, double tolerance)
{
    double violation = (B * x - z).maxCoeff();
    if(not std::isfinite(violation) or violation <= 1e-9 or violation > tolerance) return;
    for(int sweep = 0; sweep < 20 and violation > 1e-9; ++sweep)
    {
        for(Eigen::Index row = 0; row < B.rows(); ++row)
        {
            const double residual = B.row(row).dot(x) - z(row);
            const double normSquared = B.row(row).squaredNorm();
            if(residual > 1e-10 and normSquared > 1e-20)
                x -= B.row(row).transpose() * ((residual + 1e-10) / normSquared);
        }
        violation = (B * x - z).maxCoeff();
    }
}

} // namespace

int main(int argc, char **argv)
{
    if(argc != 7)
    {
        std::cerr << "Usage: e07_r0_qp_replay SNAPSHOT SEED TOL MAX_STEPS GROUP OUTPUT\n";
        return 2;
    }
    try
    {
        const std::filesystem::path directory(argv[1]);
        const std::string seedName(argv[2]);
        const double tolerance = std::stod(argv[3]);
        const unsigned int maxSteps = static_cast<unsigned int>(std::stoul(argv[4]));
        const std::string group(argv[5]);
        const std::filesystem::path output(argv[6]);
        const Eigen::MatrixXd H = read_csv(directory / "H.csv");
        const Eigen::VectorXd f = read_vector(directory / "f.csv");
        const Eigen::MatrixXd fullB = read_csv(directory / "B.csv");
        const Eigen::VectorXd fullZ = read_vector(directory / "z.csv");
        const int horizon = 12;
        std::vector<int> selected = group_rows(group, fullB.rows(), horizon);
        Eigen::MatrixXd B(selected.size(), fullB.cols());
        Eigen::VectorXd z(selected.size());
        for(size_t index = 0; index < selected.size(); ++index)
        {
            B.row(index) = fullB.row(selected[index]);
            z(index) = fullZ(selected[index]);
        }
        Eigen::VectorXd seed;
        if(seedName == "S0") seed = read_vector(directory / "U0.csv");
        else if(seedName == "S1") seed = read_vector(directory / "previous_optimal_horizon.csv");
        else if(seedName == "S2") seed = projected_seed(
            read_vector(directory / "U0.csv"), B, z);
        else if(seedName == "S3") seed = projected_seed(
            Eigen::VectorXd::Zero(H.rows()), B, z);
        else if(std::filesystem::exists(seedName)) seed = read_vector(seedName);
        else throw std::invalid_argument("Unknown seed " + seedName);
        if(seed.size() != H.rows())
            throw std::runtime_error("Seed dimension mismatch");

        SolverOptions<double> options;
        options.method = "active set";
        options.stepSizeTolerance = tolerance;
        options.maxSteps = maxSteps;
        QPSolver<double> solver(options);
        const auto start = std::chrono::steady_clock::now();
        Eigen::VectorXd solution = solver.solve(H, f, B, z, seed);
        const auto stop = std::chrono::steady_clock::now();
        const double solveTimeMs = std::chrono::duration<double, std::milli>(
            stop - start).count();
        const auto raw = solver.results();
        const double rawViolation = (B * solution - z).maxCoeff();
        Eigen::Index largestRawRow = 0;
        (B * solution - z).maxCoeff(&largestRawRow);
        const double rawObjective = 0.5 * solution.dot(H * solution) + f.dot(solution);
        write_vector(output.string() + ".raw_solution.csv", solution);
        polish(solution, B, z, tolerance);
        const double violation = (B * solution - z).maxCoeff();
        const double objective = 0.5 * solution.dot(H * solution) + f.dot(solution);
        std::vector<int> active;
        for(Eigen::Index row = 0; row < B.rows(); ++row)
            if(std::abs(B.row(row).dot(solution) - z(row)) <= 1e-7)
                active.push_back(selected[static_cast<size_t>(row)]);
        write_vector(output.string() + ".solution.csv", solution);
        std::ofstream json(output);
        json << std::setprecision(17)
             << "{\n  \"status\": \"returned\",\n"
             << "  \"seed\": \"" << seedName << "\",\n"
             << "  \"group\": \"" << group << "\",\n"
             << "  \"tolerance\": " << tolerance << ",\n"
             << "  \"max_steps\": " << maxSteps << ",\n"
             << "  \"iterations\": " << raw.numberOfSteps << ",\n"
             << "  \"final_step_size\": " << raw.finalStepSize << ",\n"
             << "  \"solve_time_ms\": " << solveTimeMs << ",\n"
             << "  \"raw_max_violation\": " << rawViolation << ",\n"
             << "  \"raw_largest_violating_row\": "
             << selected[static_cast<size_t>(largestRawRow)] << ",\n"
             << "  \"raw_objective\": " << rawObjective << ",\n"
             << "  \"polished_max_violation\": " << violation << ",\n"
             << "  \"objective\": " << objective << ",\n"
             << "  \"seed_max_violation\": " << (B * seed - z).maxCoeff() << ",\n"
             << "  \"active_constraint_rows\": [";
        for(size_t index = 0; index < active.size(); ++index)
        {
            if(index > 0) json << ',';
            json << active[index];
        }
        json << "]\n}\n";
        return 0;
    }
    catch(const std::exception &error)
    {
        std::ofstream json(argv[6]);
        json << "{\"status\":\"exception\",\"message\":\""
             << error.what() << "\"}\n";
        return 1;
    }
}
