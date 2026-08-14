/**
 * @file DiscreteIntegratorLQR.h
 * @brief Scalar discrete-LQR weights for x[k+1] = x[k] + dt u[k].
 */

#ifndef DISCRETE_INTEGRATOR_LQR_H
#define DISCRETE_INTEGRATOR_LQR_H

#include <cmath>
#include <stdexcept>

namespace RobotLibrary { namespace Math {

inline double integrator_stage_weight_for_gain(const double gain,
                                                const double controlWeight,
                                                const double dt)
{
    if(not std::isfinite(gain) or not std::isfinite(controlWeight)
       or not std::isfinite(dt) or gain <= 0.0 or controlWeight <= 0.0
       or dt <= 0.0 or gain * dt >= 1.0)
    {
        throw std::invalid_argument(
            "Discrete integrator LQR requires gain>0, R>0, dt>0 and gain*dt<1.");
    }
    return gain * gain * controlWeight / (1.0 - gain * dt);
}

inline double integrator_terminal_weight_for_gain(const double gain,
                                                   const double controlWeight,
                                                   const double dt)
{
    if(not std::isfinite(gain) or not std::isfinite(controlWeight)
       or not std::isfinite(dt) or gain <= 0.0 or controlWeight <= 0.0
       or dt <= 0.0 or gain * dt >= 1.0)
    {
        throw std::invalid_argument(
            "Discrete integrator LQR requires gain>0, R>0, dt>0 and gain*dt<1.");
    }
    return gain * controlWeight / (dt * (1.0 - gain * dt));
}

inline double integrator_terminal_weight(const double stateWeight,
                                          const double controlWeight,
                                          const double dt)
{
    if(not std::isfinite(stateWeight) or not std::isfinite(controlWeight)
       or not std::isfinite(dt) or stateWeight <= 0.0
       or controlWeight <= 0.0 or dt <= 0.0)
    {
        throw std::invalid_argument(
            "Discrete integrator LQR requires Q>0, R>0 and dt>0.");
    }
    return 0.5 * (stateWeight + std::sqrt(
        stateWeight * stateWeight
        + 4.0 * stateWeight * controlWeight / (dt * dt)));
}

} } // namespace RobotLibrary::Math

#endif
