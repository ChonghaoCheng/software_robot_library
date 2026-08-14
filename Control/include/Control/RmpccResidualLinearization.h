/**
 * @file RmpccResidualLinearization.h
 * @brief Complete finite-difference linearisation of FullScrewSE3 residuals.
 */

#ifndef RMPCC_RESIDUAL_LINEARIZATION_H
#define RMPCC_RESIDUAL_LINEARIZATION_H

#include <Control/RmpccCostGeometry.h>
#include <Control/RmpccPrediction.h>

#include <Eigen/Core>

#include <cmath>
#include <stdexcept>

namespace RobotLibrary { namespace Control {

enum class RmpccResidualLinearization
{
    FrozenProjector,
    FullResidualJacobian
};

struct RmpccFullScrewResiduals
{
    Eigen::Vector<double,6> contour = Eigen::Vector<double,6>::Zero();
    Eigen::Vector<double,6> lag = Eigen::Vector<double,6>::Zero();
};

struct RmpccFullScrewResidualLinearization
{
    RmpccFullScrewResiduals residual;
    Eigen::Matrix<double,6,7> contourJacobian =
        Eigen::Matrix<double,6,7>::Zero();
    Eigen::Matrix<double,6,7> lagJacobian =
        Eigen::Matrix<double,6,7>::Zero();
};

template<typename ReferenceTangentFunction>
RmpccFullScrewResiduals
rmpcc_full_screw_residuals(
    const RmpccStateVector &state,
    const Eigen::Matrix<double,6,6> &metric,
    const double regularization,
    ReferenceTangentFunction &&referenceTangent)
{
    const Eigen::Vector<double,6> tangent = referenceTangent(state(6));
    const Eigen::Vector<double,6> errorTangent =
        rmpcc_error_coordinate_path_tangent(state.head<6>(), tangent);
    const RmpccErrorProjection projection = rmpcc_error_projection(
        errorTangent, metric, RmpccLagGeometry::FullScrew, regularization);

    RmpccFullScrewResiduals result;
    result.contour = projection.contour * state.head<6>();
    result.lag = projection.lag * state.head<6>();
    return result;
}

template<typename ReferenceTangentFunction>
RmpccFullScrewResidualLinearization
rmpcc_linearize_full_screw_residuals(
    const RmpccStateVector &state,
    const Eigen::Matrix<double,6,6> &metric,
    const double regularization,
    const double finiteDifferenceStep,
    ReferenceTangentFunction &&referenceTangent)
{
    if(not std::isfinite(finiteDifferenceStep) or finiteDifferenceStep <= 0.0)
    {
        throw std::invalid_argument(
            "[ERROR] [RMPCC RESIDUAL] finiteDifferenceStep must be positive.");
    }

    auto &&tangent = referenceTangent;
    RmpccFullScrewResidualLinearization result;
    result.residual = rmpcc_full_screw_residuals(
        state, metric, regularization, tangent);
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
        if(std::abs(denominator) <= 1e-15)
        {
            continue;
        }
        const RmpccFullScrewResiduals plusResidual =
            rmpcc_full_screw_residuals(plus, metric, regularization, tangent);
        const RmpccFullScrewResiduals minusResidual =
            rmpcc_full_screw_residuals(minus, metric, regularization, tangent);
        result.contourJacobian.col(column) =
            (plusResidual.contour - minusResidual.contour)
            / denominator;
        result.lagJacobian.col(column) =
            (plusResidual.lag - minusResidual.lag)
            / denominator;
    }
    return result;
}

} } // namespace RobotLibrary::Control

#endif
