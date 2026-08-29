#ifndef ROBOT_LIBRARY_MATH_QP_ACCEPTANCE_H
#define ROBOT_LIBRARY_MATH_QP_ACCEPTANCE_H

#include <Math/QPSolver.h>

#include <cmath>
#include <stdexcept>
#include <string>

namespace RobotLibrary { namespace Math {

inline const char *solver_termination_reason_name(const SolverTerminationReason reason)
{
    return reason == SolverTerminationReason::Converged
        ? "Converged" : "MaxIterations";
}

/** Require explicit convergence before a QP solution may become a command. */
template<typename DataType>
inline void require_qp_converged(const SolverResults<DataType> &results,
                                 const std::string &context)
{
    if(results.terminationReason != SolverTerminationReason::Converged)
    {
        throw std::runtime_error(
            context + ": QP did not converge (termination="
            + solver_termination_reason_name(results.terminationReason)
            + ", iterations=" + std::to_string(results.numberOfSteps)
            + ", final_step=" + std::to_string(results.finalStepSize) + ").");
    }
}

template<typename DataType, typename Derived>
inline void require_qp_result_accepted(
    const SolverResults<DataType> &results,
    const Eigen::MatrixBase<Derived> &solution,
    const DataType primalViolation,
    const DataType primalTolerance,
    const std::string &context)
{
    require_qp_converged(results, context);
    if(solution.size() == 0 || !solution.allFinite())
        throw std::runtime_error(context + ": QP returned a non-finite/empty solution.");
    if(!std::isfinite(primalViolation) || primalViolation > primalTolerance)
        throw std::runtime_error(context + ": QP returned an infeasible solution.");
}

}} // namespace RobotLibrary::Math

#endif
