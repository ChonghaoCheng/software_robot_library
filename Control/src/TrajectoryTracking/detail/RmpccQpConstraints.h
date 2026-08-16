/**
 * @file TrajectoryTracking/detail/RmpccQpConstraints.h
 * @brief Pure construction helpers for RMPCC QP constraints.
 */

#ifndef RMPCC_QP_CONSTRAINTS_H
#define RMPCC_QP_CONSTRAINTS_H

#include "RmpccProgressConstraints.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace RobotLibrary { namespace Control {

struct RmpccQpConstraints
{
    Eigen::MatrixXd Aeq;
    Eigen::VectorXd yeq;
    Eigen::MatrixXd Bineq;
    Eigen::VectorXd zineq;
};

inline bool rmpcc_fixed_progress_schedule_satisfies_limits(
    const Eigen::VectorXd &fixedProgressRates,
    const double dt,
    const double completionRemaining,
    const double scheduleRemaining,
    const double tolerance = 1e-12)
{
    if(fixedProgressRates.size() == 0)
    {
        return true;
    }
    if(not fixedProgressRates.allFinite() or not std::isfinite(dt)
       or dt <= 0.0 or not std::isfinite(completionRemaining)
       or not std::isfinite(scheduleRemaining))
    {
        return false;
    }
    const double completionAdvance = dt * fixedProgressRates.sum();
    const double scheduledAdvance = dt * fixedProgressRates(0);
    return completionAdvance <= completionRemaining + tolerance
        and scheduledAdvance <= scheduleRemaining + tolerance;
}

inline RmpccQpConstraints rmpcc_build_qp_constraints(
    const Eigen::VectorXd &lower,
    const Eigen::VectorXd &upper,
    const int horizon,
    const double dt,
    const double completionRemaining,
    const double scheduleRemaining,
    const bool fixedProgressSchedule,
    const Eigen::VectorXd &fixedProgressRates)
{
    if(horizon < 1 or lower.size() != upper.size()
       or lower.size() != 7 * horizon
       or fixedProgressRates.size() != horizon)
    {
        throw std::invalid_argument(
            "RMPCC QP constraint dimensions do not match the horizon.");
    }

    const Eigen::Index variableDim = lower.size();
    const Eigen::Index progressOffset = variableDim - horizon;
    RmpccQpConstraints constraints;
    if(fixedProgressSchedule)
    {
        if(not rmpcc_fixed_progress_schedule_satisfies_limits(
               fixedProgressRates, dt, completionRemaining, scheduleRemaining))
        {
            throw std::invalid_argument(
                "Fixed RMPCC progress schedule violates completion or supplied-schedule limits.");
        }
        constraints.Aeq = Eigen::MatrixXd::Zero(horizon, variableDim);
        constraints.yeq = fixedProgressRates;
        for(int stage = 0; stage < horizon; ++stage)
        {
            constraints.Aeq(stage, progressOffset + stage) = 1.0;
        }

        // Fixed progress is represented only by independent equalities. The
        // completion and first-step schedule rows are algebraically implied by
        // those equalities and were checked above.
        constraints.Bineq = Eigen::MatrixXd::Zero(2 * progressOffset, variableDim);
        constraints.zineq = Eigen::VectorXd::Zero(2 * progressOffset);
        constraints.Bineq.topLeftCorner(progressOffset, progressOffset).setIdentity();
        constraints.zineq.head(progressOffset) = upper.head(progressOffset);
        constraints.Bineq.bottomLeftCorner(progressOffset, progressOffset) =
            -Eigen::MatrixXd::Identity(progressOffset, progressOffset);
        constraints.zineq.tail(progressOffset) = -lower.head(progressOffset);
    }
    else
    {
        constraints.Aeq = Eigen::MatrixXd(0, variableDim);
        constraints.yeq = Eigen::VectorXd(0);
        constraints.Bineq = Eigen::MatrixXd::Zero(
            2 * variableDim + 2, variableDim);
        constraints.zineq = Eigen::VectorXd::Zero(2 * variableDim + 2);
        constraints.Bineq.topRows(variableDim).setIdentity();
        constraints.zineq.head(variableDim) = upper;
        constraints.Bineq.block(variableDim, 0, variableDim, variableDim) =
            -Eigen::MatrixXd::Identity(variableDim, variableDim);
        constraints.zineq.segment(variableDim, variableDim) = -lower;
        constraints.Bineq.block(2 * variableDim, progressOffset, 1, horizon) =
            rmpcc_completion_progress_row(horizon, dt).transpose();
        constraints.Bineq.block(2 * variableDim + 1, progressOffset, 1, horizon) =
            rmpcc_schedule_progress_row(horizon, dt).transpose();
        constraints.zineq(2 * variableDim) = completionRemaining;
        constraints.zineq(2 * variableDim + 1) = scheduleRemaining;
    }
    return constraints;
}

} } // namespace RobotLibrary::Control

#endif
