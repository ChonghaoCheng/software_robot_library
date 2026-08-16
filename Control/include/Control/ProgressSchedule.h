#ifndef CONTROL_PROGRESS_SCHEDULE_H
#define CONTROL_PROGRESS_SCHEDULE_H

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace RobotLibrary { namespace Control {

inline double reference_schedule_progress(
    const std::uint64_t controlStep,
    const double controlDt,
    const double trajectoryDuration)
{
    if(not std::isfinite(controlDt) or controlDt <= 0.0
       or not std::isfinite(trajectoryDuration) or trajectoryDuration <= 0.0)
    {
        throw std::invalid_argument("reference schedule requires positive finite dt and duration");
    }
    return std::clamp(
        static_cast<double>(controlStep) * controlDt / trajectoryDuration,
        0.0, 1.0);
}

inline Eigen::VectorXd reference_schedule_rates(
    const double initialProgress,
    const int horizonSteps,
    const double controlDt,
    const double trajectoryDuration)
{
    if(horizonSteps < 0 or not std::isfinite(initialProgress)
       or not std::isfinite(controlDt) or controlDt <= 0.0
       or not std::isfinite(trajectoryDuration) or trajectoryDuration <= 0.0)
    {
        throw std::invalid_argument("reference schedule requires valid progress, horizon, dt and duration");
    }
    Eigen::VectorXd rates(horizonSteps);
    double progress = std::clamp(initialProgress, 0.0, 1.0);
    const double nominalRate = 1.0 / trajectoryDuration;
    for(int stage = 0; stage < horizonSteps; ++stage)
    {
        const double next = std::clamp(progress + controlDt * nominalRate, 0.0, 1.0);
        rates(stage) = (next - progress) / controlDt;
        progress = next;
    }
    return rates;
}

} } // namespace RobotLibrary::Control

#endif
