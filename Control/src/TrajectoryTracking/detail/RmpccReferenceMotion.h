/**
 * @file RmpccReferenceMotion.h
 * @brief Experimental separation of geometric path tangent and finite stage motion.
 */

#ifndef RMPCC_REFERENCE_MOTION_H
#define RMPCC_REFERENCE_MOTION_H

#include "RmpccPrediction.h"

#include <Math/MathFunctions.h>

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace RobotLibrary { namespace Control {

struct RmpccPathVelocityResidualLinearization
{
    Eigen::Vector<double,6> residual = Eigen::Vector<double,6>::Zero();
    Eigen::Matrix<double,6,7> stateJacobian = Eigen::Matrix<double,6,7>::Zero();
    Eigen::Matrix<double,6,7> inputJacobian = Eigen::Matrix<double,6,7>::Zero();
};

/** Centred body-frame Lie tangent dT/ds, expressed at T(s). */
template<typename ReferenceTransformFunction>
Eigen::Vector<double,6>
rmpcc_centred_geometric_tangent(
    const double progress,
    const double epsilon,
    ReferenceTransformFunction &&referenceTransform)
{
    if(not std::isfinite(progress) or not std::isfinite(epsilon) or epsilon <= 0.0)
    {
        throw std::invalid_argument(
            "[ERROR] [RMPCC REFERENCE MOTION] finite progress and positive epsilon required.");
    }
    auto &&reference = referenceTransform;
    const double s = std::clamp(progress, 0.0, 1.0);
    const double sMinus = std::clamp(s - epsilon, 0.0, 1.0);
    const double sPlus = std::clamp(s + epsilon, 0.0, 1.0);
    const double denominator = sPlus - sMinus;
    if(denominator <= 1e-15)
    {
        return Eigen::Vector<double,6>::Zero();
    }
    const Eigen::Matrix4d current = reference(s);
    const Eigen::Matrix4d currentInverse = RobotLibrary::Math::se3_inverse(current);
    const Eigen::Vector<double,6> plus = RobotLibrary::Math::se3_logarithm(
        currentInverse * reference(sPlus));
    const Eigen::Vector<double,6> minus = RobotLibrary::Math::se3_logarithm(
        currentInverse * reference(sMinus));
    return (plus - minus) / denominator;
}

/** Finite body-frame reference motion over one stage, expressed per unit time. */
template<typename ReferenceTransformFunction>
Eigen::Vector<double,6>
rmpcc_stage_reference_motion(
    const double progress,
    const double progressRate,
    const double dt,
    ReferenceTransformFunction &&referenceTransform)
{
    if(not std::isfinite(progress) or not std::isfinite(progressRate)
       or not std::isfinite(dt) or dt <= 0.0)
    {
        throw std::invalid_argument(
            "[ERROR] [RMPCC REFERENCE MOTION] finite state/rate and positive dt required.");
    }
    auto &&reference = referenceTransform;
    const double s = std::clamp(progress, 0.0, 1.0);
    const double next = std::clamp(s + dt * progressRate, 0.0, 1.0);
    return RobotLibrary::Math::se3_logarithm(
               RobotLibrary::Math::se3_inverse(reference(s)) * reference(next)) / dt;
}

/** Error-frame feedforward that exactly realizes the selected finite stage motion. */
template<typename ReferenceTransformFunction>
Eigen::Vector<double,6>
rmpcc_stage_consistent_feedforward(
    const Eigen::Vector<double,6> &error,
    const double progress,
    const double progressRate,
    const double dt,
    ReferenceTransformFunction &&referenceTransform)
{
    const Eigen::Matrix4d relative = RobotLibrary::Math::se3_exponential(error);
    return RobotLibrary::Math::adjoint(RobotLibrary::Math::se3_inverse(relative))
           * rmpcc_stage_reference_motion(
               progress, progressRate, dt,
               std::forward<ReferenceTransformFunction>(referenceTransform));
}

template<typename ReferenceTransformFunction>
Eigen::Vector<double,6>
rmpcc_stage_consistent_path_velocity_residual(
    const RmpccStateVector &state,
    const RmpccInputVector &input,
    const double dt,
    ReferenceTransformFunction &&referenceTransform)
{
    return input.head<6>() - rmpcc_stage_consistent_feedforward(
        state.head<6>(), state(6), input(6), dt,
        std::forward<ReferenceTransformFunction>(referenceTransform));
}

/**
 * Complete local linearization of r(E,s,u,v)=u-Ad(E^-1)xi_step(s,v).
 * J_u is exact identity. Error, progress and progress-rate derivatives use
 * bounded centred finite differences at the RTI nominal point.
 */
template<typename ReferenceTransformFunction>
RmpccPathVelocityResidualLinearization
rmpcc_linearize_stage_consistent_path_velocity_residual(
    const RmpccStateVector &state,
    const RmpccInputVector &input,
    const double dt,
    const double finiteDifferenceStep,
    ReferenceTransformFunction &&referenceTransform)
{
    if(not std::isfinite(finiteDifferenceStep) or finiteDifferenceStep <= 0.0)
    {
        throw std::invalid_argument(
            "[ERROR] [RMPCC REFERENCE MOTION] finiteDifferenceStep must be positive.");
    }
    auto &&reference = referenceTransform;
    RmpccPathVelocityResidualLinearization result;
    result.residual = rmpcc_stage_consistent_path_velocity_residual(
        state, input, dt, reference);

    for(int column = 0; column < 7; ++column)
    {
        RmpccStateVector plus = state;
        RmpccStateVector minus = state;
        plus(column) += finiteDifferenceStep;
        minus(column) -= finiteDifferenceStep;
        double denominator = 2.0 * finiteDifferenceStep;
        if(column == 6)
        {
            plus(6) = std::clamp(plus(6), 0.0, 1.0);
            minus(6) = std::clamp(minus(6), 0.0, 1.0);
            denominator = plus(6) - minus(6);
        }
        if(std::abs(denominator) > 1e-15)
        {
            result.stateJacobian.col(column) =
                (rmpcc_stage_consistent_path_velocity_residual(
                     plus, input, dt, reference)
                 - rmpcc_stage_consistent_path_velocity_residual(
                     minus, input, dt, reference)) / denominator;
        }
    }

    result.inputJacobian.leftCols<6>().setIdentity();
    const auto ratePerturbations = rmpcc_progress_rate_perturbations(
        state(6), input(6), dt, finiteDifferenceStep);
    const double rateDenominator = ratePerturbations.second - ratePerturbations.first;
    if(std::abs(rateDenominator) > 1e-15)
    {
        RmpccInputVector plus = input;
        RmpccInputVector minus = input;
        plus(6) = ratePerturbations.second;
        minus(6) = ratePerturbations.first;
        result.inputJacobian.col(6) =
            (rmpcc_stage_consistent_path_velocity_residual(
                 state, plus, dt, reference)
             - rmpcc_stage_consistent_path_velocity_residual(
                 state, minus, dt, reference)) / rateDenominator;
    }
    return result;
}

} } // namespace RobotLibrary::Control

#endif
