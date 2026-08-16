/**
 * @file TrajectoryTracking/detail/RmpccProgressConstraints.h
 * @brief Pure progress-constraint helpers shared by RMPCC and its tests.
 */

#ifndef RMPCC_PROGRESS_CONSTRAINTS_H
#define RMPCC_PROGRESS_CONSTRAINTS_H

#include <Eigen/Core>

#include <algorithm>

namespace RobotLibrary { namespace Control {

inline Eigen::VectorXd
rmpcc_completion_progress_row(const int horizon, const double dt)
{
    return Eigen::VectorXd::Constant(horizon, dt);
}

inline Eigen::VectorXd
rmpcc_schedule_progress_row(const int horizon, const double dt)
{
    Eigen::VectorXd row = Eigen::VectorXd::Zero(horizon);
    if(horizon > 0)
    {
        row(0) = dt;
    }
    return row;
}

inline Eigen::VectorXd
rmpcc_progress_rate_lower_bounds(const int horizon,
                                 const double progressRateMin,
                                 const bool relaxForCompletion,
                                 const double dt,
                                 const double scheduleRemaining)
{
    Eigen::VectorXd lower = Eigen::VectorXd::Constant(
        horizon, relaxForCompletion ? 0.0 : progressRateMin);
    if(not relaxForCompletion and horizon > 0
       and scheduleRemaining < dt * progressRateMin)
    {
        lower(0) = 0.0;
    }
    return lower;
}

inline void
rmpcc_clip_progress_rates(Eigen::VectorXd &rates,
                          const Eigen::VectorXd &lower,
                          const Eigen::VectorXd &upper,
                          const double dt,
                          const double completionRemaining,
                          const double scheduleRemaining)
{
    for(int stage = 0; stage < rates.size(); ++stage)
    {
        rates(stage) = std::clamp(rates(stage), lower(stage), upper(stage));
    }

    const double totalAdvance = dt * rates.sum();
    if(totalAdvance > completionRemaining && rates.size() > 0)
    {
        const double feasibleRate = std::max(0.0, completionRemaining)
            / (static_cast<double>(rates.size()) * std::max(dt, 1e-9));
        for(int stage = 0; stage < rates.size(); ++stage)
        {
            rates(stage) = std::clamp(feasibleRate, lower(stage), upper(stage));
        }
    }

    if(rates.size() > 0)
    {
        rates(0) = std::min(
            rates(0), std::max(0.0, scheduleRemaining) / std::max(dt, 1e-9));
    }
}

} } // namespace RobotLibrary::Control

#endif
