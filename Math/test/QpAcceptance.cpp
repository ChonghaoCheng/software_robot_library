#include <Math/QpAcceptance.h>

#include <Eigen/Core>

#include <iostream>
#include <limits>
#include <stdexcept>

using RobotLibrary::Math::require_qp_result_accepted;

namespace {
bool rejects(const SolverResults<double> &results,
             const Eigen::VectorXd &solution,
             const double violation)
{
    try
    {
        require_qp_result_accepted(results, solution, violation, 1e-6, "test");
        return false;
    }
    catch(const std::runtime_error &)
    {
        return true;
    }
}
}

int main()
{
    SolverResults<double> result;
    result.terminationReason = SolverTerminationReason::Converged;
    result.numberOfSteps = 3;
    result.finalStepSize = 1e-8;
    const Eigen::VectorXd finite = Eigen::Vector2d(0.1, -0.2);
    require_qp_result_accepted(result, finite, 0.0, 1e-6, "converged");

    result.terminationReason = SolverTerminationReason::MaxIterations;
    if(!rejects(result, finite, 0.0))
        throw std::runtime_error("feasible MaxIterations result was accepted");
    if(!rejects(result, finite, -1e-9))
        throw std::runtime_error("strictly feasible MaxIterations result was accepted");

    result.terminationReason = SolverTerminationReason::Converged;
    Eigen::VectorXd nonfinite = finite;
    nonfinite(0) = std::numeric_limits<double>::quiet_NaN();
    if(!rejects(result, nonfinite, 0.0))
        throw std::runtime_error("non-finite result was accepted");
    if(!rejects(result, finite, 1e-3))
        throw std::runtime_error("infeasible result was accepted");

    result.terminationReason = SolverTerminationReason::MaxIterations;
    try
    {
        require_qp_result_accepted(result, finite, 0.0, 1e-6, "preserve-reason");
    }
    catch(const std::runtime_error &error)
    {
        if(std::string(error.what()).find("MaxIterations") == std::string::npos)
            throw std::runtime_error("termination reason was not preserved");
        std::cout << "qp_acceptance_test PASS\n";
        return 0;
    }
    throw std::runtime_error("MaxIterations did not throw");
}
