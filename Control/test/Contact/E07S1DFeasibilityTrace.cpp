/** Offline E07-S1D trace harness for the real RobotLibrary active-set path. */

#include <Math/QPSolver.h>

#include <Eigen/Core>

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
    Eigen::MatrixXd result(rows.size(), rows.empty() ? 0 : rows.front().size());
    for(Eigen::Index row = 0; row < result.rows(); ++row)
        for(Eigen::Index column = 0; column < result.cols(); ++column)
            result(row, column) = rows[static_cast<size_t>(row)]
                                      [static_cast<size_t>(column)];
    return result;
}

Eigen::VectorXd read_vector(const std::filesystem::path &path)
{
    const Eigen::MatrixXd data = read_csv(path);
    if(data.cols() != 1) throw std::runtime_error("Expected vector " + path.string());
    return data.col(0);
}

template<class Derived>
void write_semicolon_vector(std::ostream &stream,
                            const Eigen::MatrixBase<Derived> &values)
{
    for(Eigen::Index index = 0; index < values.size(); ++index)
    {
        if(index) stream << ';';
        stream << values(index);
    }
}

void write_int_vector(std::ostream &stream, const std::vector<int> &values)
{
    for(size_t index = 0; index < values.size(); ++index)
    {
        if(index) stream << ';';
        stream << values[index];
    }
}

void write_vector(const std::filesystem::path &path, const Eigen::VectorXd &x)
{
    std::ofstream stream(path);
    stream << std::setprecision(17);
    for(Eigen::Index index = 0; index < x.size(); ++index) stream << x(index) << '\n';
}

} // namespace

int main(int argc, char **argv)
{
    if(argc != 4)
    {
        std::cerr << "usage: e07_s1d_feasibility_trace SNAPSHOT SEED OUTPUT\n";
        return 2;
    }
    try
    {
        const std::filesystem::path snapshot(argv[1]);
        const std::filesystem::path seedPath(argv[2]);
        const std::filesystem::path output(argv[3]);
        std::filesystem::create_directories(output);
        const Eigen::MatrixXd H = read_csv(snapshot / "H.csv");
        const Eigen::VectorXd f = read_vector(snapshot / "f.csv");
        const Eigen::MatrixXd B = read_csv(snapshot / "B.csv");
        const Eigen::VectorXd z = read_vector(snapshot / "z.csv");
        const Eigen::VectorXd seed = read_vector(seedPath);

        SolverOptions<double> options;
        options.method = "active set";
        options.stepSizeTolerance = 1e-4;
        options.maxSteps = 50;
        QPSolver<double> solver(options);

        std::ofstream trace(output / "raw_trace.csv");
        trace << std::setprecision(17)
              << "sequence,event,iteration,step_size,alpha,blocking_constraint,"
                 "ratio_constraint,ratio_slack,ratio_added_distance,raw_ratio,"
                 "schur_solve_residual,schur_rhs_norm,schur_acceptance_threshold,"
                 "schur_ldlt_accepted,schur_cod_fallback,"
                 "active_set,inactive_set,x,dx\n";
        size_t sequence = 0;
        solver.set_active_set_trace_sink(
            [&](const ActiveSetTraceRecord<double> &record)
            {
                trace << sequence++ << ',' << record.event << ',' << record.iteration
                      << ',' << record.stepSize << ',' << record.alpha
                      << ',' << record.blockingConstraint
                      << ',' << record.ratioConstraint
                      << ',' << record.ratioSlack
                      << ',' << record.ratioAddedDistance
                      << ',' << record.rawRatio
                      << ',' << record.schurSolveResidual
                      << ',' << record.schurRhsNorm
                      << ',' << record.schurAcceptanceThreshold
                      << ',' << (record.schurLdltAccepted ? 1 : 0)
                      << ',' << (record.schurCodFallbackUsed ? 1 : 0) << ',';
                write_int_vector(trace, record.activeSet); trace << ',';
                write_int_vector(trace, record.inactiveSet); trace << ',';
                write_semicolon_vector(trace, record.x); trace << ',';
                write_semicolon_vector(trace, record.dx); trace << '\n';
            });
        const Eigen::VectorXd solution = solver.solve(H, f, B, z, seed);
        trace.close();
        write_vector(output / "raw_solution.csv", solution);
        const auto result = solver.results();
        const Eigen::VectorXd residual = B * solution - z;
        Eigen::Index largest = 0;
        const double rp = residual.maxCoeff(&largest);
        const double objective = 0.5 * solution.dot(H * solution) + f.dot(solution);
        std::ofstream summary(output / "raw_summary.json");
        summary << std::setprecision(17)
                << "{\n"
                << "  \"converged\": " << (result.converged ? "true" : "false") << ",\n"
                << "  \"iterations\": " << result.numberOfSteps << ",\n"
                << "  \"max_steps\": 50,\n"
                << "  \"step_size_tolerance\": 0.0001,\n"
                << "  \"final_step_size\": " << result.finalStepSize << ",\n"
                << "  \"seed_rp\": " << (B * seed - z).maxCoeff() << ",\n"
                << "  \"raw_rp\": " << rp << ",\n"
                << "  \"largest_violating_row\": " << largest << ",\n"
                << "  \"objective\": " << objective << "\n"
                << "}\n";
        return 0;
    }
    catch(const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
