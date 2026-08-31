#include <Control/Core/DataStructures.h>
#include <Control/TrajectoryTracking/SerialLinkRMPCC.h>
#include <Math/BoxAwareActiveSet.h>

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <string>

namespace {

bool check(const bool condition, const std::string &message)
{
    if(not condition)
    {
        std::cerr << "FAIL: " << message << '\n';
    }
    return condition;
}

bool budget_isolation_contract()
{
    RobotLibrary::Control::SerialLinkParameters parameters;
    parameters.qpsolver.maxSteps = 50;
    const bool fallback = check(
        parameters.optimized_progress_qp_max_steps() == 50,
        "absent task budget falls back to shared budget");
    parameters.taskQpMaxSteps = 75;
    return fallback
        and check(parameters.optimized_progress_qp_max_steps() == 75,
                  "task budget resolves independently")
        and check(parameters.qpsolver.maxSteps == 50,
                  "resolved-rate budget remains unchanged");
}

bool direct_objective_telemetry_contract()
{
    Eigen::Matrix2d H = Eigen::Matrix2d::Identity();
    Eigen::Vector2d f;
    f << -0.25, 0.5;
    Eigen::Matrix<double, 4, 2> B;
    B << 1.0, 0.0,
         0.0, 1.0,
        -1.0, 0.0,
         0.0,-1.0;
    Eigen::Vector4d z = Eigen::Vector4d::Ones();
    Eigen::Vector2d seed = Eigen::Vector2d::Zero();
    SolverOptions<double> options;
    options.maxSteps = 75;
    options.stepSizeTolerance = 1e-5;

    const auto result = solve_box_aware_active_set(H, f, B, z, seed, options);
    RobotLibrary::Control::RmpccDiagnostics diagnostics;
    diagnostics.record_qp_objective(result.solver);

    return check(result.solver.terminationReason == SolverTerminationReason::Converged,
                 "representative optimized-progress QP converges")
        and check(std::isfinite(result.solver.objectiveFunction),
                  "solver objective is finite")
        and check(diagnostics.qpObjective == result.solver.objectiveFunction,
                  "diagnostic is assigned directly from solver result");
}

} // namespace

int main()
{
    const bool passed = budget_isolation_contract()
        and direct_objective_telemetry_contract();
    if(passed)
    {
        std::cout << "PASS: task-QP budget isolation and RMPCC objective telemetry\n";
        return 0;
    }
    return 1;
}
