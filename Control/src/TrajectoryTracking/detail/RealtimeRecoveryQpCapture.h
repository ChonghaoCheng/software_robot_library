#ifndef ROBOT_LIBRARY_REALTIME_RECOVERY_QP_CAPTURE_H
#define ROBOT_LIBRARY_REALTIME_RECOVERY_QP_CAPTURE_H

#include <Math/QPSolver.h>

#include <Eigen/Dense>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <string>

namespace RobotLibrary { namespace Control {

inline void realtime_recovery_write_matrix(
    const std::filesystem::path &path, const Eigen::Ref<const Eigen::MatrixXd> &matrix)
{
    std::ofstream stream(path);
    if(!stream) throw std::runtime_error("Unable to write frozen QP matrix: " + path.string());
    stream << std::setprecision(17);
    for(int row = 0; row < matrix.rows(); ++row)
    {
        for(int column = 0; column < matrix.cols(); ++column)
        {
            if(column) stream << ',';
            stream << matrix(row, column);
        }
        stream << '\n';
    }
}

inline void realtime_recovery_write_vector(
    const std::filesystem::path &path, const Eigen::Ref<const Eigen::VectorXd> &vector)
{
    realtime_recovery_write_matrix(path, vector);
}

/** Opt-in, synchronous frozen-QP capture. Never enabled in timing or production runs. */
inline void write_realtime_recovery_qp_capture(
    const std::string &controller,
    const std::uint64_t controlStep,
    const Eigen::Ref<const Eigen::MatrixXd> &H,
    const Eigen::Ref<const Eigen::VectorXd> &f,
    const Eigen::Ref<const Eigen::MatrixXd> &A,
    const Eigen::Ref<const Eigen::VectorXd> &y,
    const Eigen::Ref<const Eigen::MatrixXd> &B,
    const Eigen::Ref<const Eigen::VectorXd> &z,
    const Eigen::Ref<const Eigen::VectorXd> &seed,
    const Eigen::Ref<const Eigen::VectorXd> &solution,
    const SolverResults<double> &results)
{
    const char *rootText = std::getenv("ROBOT_LIBRARY_REALTIME_RECOVERY_QP_CAPTURE_DIR");
    if(rootText == nullptr || *rootText == '\0') return;

    static std::map<std::string, unsigned int> captures;
    const bool regular = controlStep == 1000 || controlStep == 5000 || controlStep == 10000;
    const bool iterationLimited =
        results.terminationReason == SolverTerminationReason::MaxIterations;
    if((!regular && !iterationLimited) || captures[controller] >= 16) return;

    const std::string reason = iterationLimited ? "max_iterations" : "regular";
    const std::filesystem::path directory = std::filesystem::path(rootText) / controller /
        ("step_" + std::to_string(controlStep) + "__" + reason);
    std::filesystem::create_directories(directory);
    realtime_recovery_write_matrix(directory / "H.csv", H);
    realtime_recovery_write_vector(directory / "f.csv", f);
    realtime_recovery_write_matrix(directory / "A.csv", A);
    realtime_recovery_write_vector(directory / "y.csv", y);
    realtime_recovery_write_matrix(directory / "B.csv", B);
    realtime_recovery_write_vector(directory / "z.csv", z);
    realtime_recovery_write_vector(directory / "seed.csv", seed);
    realtime_recovery_write_vector(directory / "production_solution.csv", solution);
    std::ofstream metadata(directory / "metadata.json");
    metadata << std::setprecision(17)
             << "{\n"
             << "  \"schema\": \"realtime_recovery_02_frozen_qp_v1\",\n"
             << "  \"controller\": \"" << controller << "\",\n"
             << "  \"control_step_index_zero_based\": " << controlStep << ",\n"
             << "  \"production_number_of_steps\": " << results.numberOfSteps << ",\n"
             << "  \"production_termination_reason\": \""
             << (iterationLimited ? "MaxIterations" : "Converged") << "\",\n"
             << "  \"production_final_step_size\": " << results.finalStepSize << ",\n"
             << "  \"production_objective\": " << results.objectiveFunction << ",\n"
             << "  \"active_set_changes\": " << results.activeSetChanges << ",\n"
             << "  \"maximum_active_set_size\": " << results.maximumActiveSetSize << ",\n"
             << "  \"unique_active_constraints\": " << results.uniqueActiveConstraints << "\n"
             << "}\n";
    ++captures[controller];
}

}} // namespace RobotLibrary::Control

#endif
