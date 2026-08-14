/**
 * @file RmpccCostGeometry.h
 * @brief Pure contour/lag projection geometry used by RMPCC running costs.
 */

#ifndef RMPCC_COST_GEOMETRY_H
#define RMPCC_COST_GEOMETRY_H

#include <Math/MathFunctions.h>
#include <Eigen/Core>

#include <algorithm>
#include <cmath>

namespace RobotLibrary { namespace Control {

enum class RmpccLagGeometry
{
    FullScrew,
    SplitTranslationRotation,
    TranslationOnly,
    RotationOnly
};

enum class RmpccObjectiveGeometry
{
    FullScrewSE3,
    DecoupledCartesianSO3
};

struct RmpccErrorProjection
{
    Eigen::Matrix<double,6,6> lag = Eigen::Matrix<double,6,6>::Zero();
    Eigen::Matrix<double,6,6> contour = Eigen::Matrix<double,6,6>::Identity();
};

/** Scale translation and rotation blocks of a six-dimensional quadratic cost. */
inline Eigen::Matrix<double,6,6>
rmpcc_component_scaled_weight(
    const Eigen::Matrix<double,6,6> &weight,
    const double translationScale,
    const double rotationScale)
{
    Eigen::Matrix<double,6,6> scale = Eigen::Matrix<double,6,6>::Identity();
    scale.diagonal().head<3>().setConstant(std::sqrt(translationScale));
    scale.diagonal().tail<3>().setConstant(std::sqrt(rotationScale));
    return scale * weight * scale;
}

template<int Dimension>
inline Eigen::Matrix<double,Dimension,Dimension>
rmpcc_metric_projection(
    const Eigen::Matrix<double,Dimension,1> &tangent,
    const Eigen::Matrix<double,Dimension,Dimension> &metric,
    const double regularization)
{
    const double denominator =
        std::max((tangent.transpose() * metric * tangent)(0), regularization);
    return tangent * (tangent.transpose() * metric) / denominator;
}

/**
 * @brief Construct the lag projector used at a non-terminal RMPCC stage.
 *
 * FullScrew reproduces the original single scalar screw-phase projection.
 * The split modes remove its translation/rotation cross blocks. A zero path
 * tangent in either subspace produces a zero lag projector there.
 */
inline RmpccErrorProjection
rmpcc_error_projection(
    const Eigen::Matrix<double,6,1> &errorTangent,
    const Eigen::Matrix<double,6,6> &metric,
    const RmpccLagGeometry geometry,
    const double regularization)
{
    RmpccErrorProjection result;
    if(geometry == RmpccLagGeometry::FullScrew)
    {
        result.lag = rmpcc_metric_projection<6>(
            errorTangent, metric, regularization);
    }
    else
    {
        if(geometry == RmpccLagGeometry::SplitTranslationRotation
           or geometry == RmpccLagGeometry::TranslationOnly)
        {
            result.lag.block<3,3>(0,0) = rmpcc_metric_projection<3>(
                errorTangent.head<3>(), metric.block<3,3>(0,0), regularization);
        }
        if(geometry == RmpccLagGeometry::SplitTranslationRotation
           or geometry == RmpccLagGeometry::RotationOnly)
        {
            result.lag.block<3,3>(3,3) = rmpcc_metric_projection<3>(
                errorTangent.tail<3>(), metric.block<3,3>(3,3), regularization);
        }
    }
    result.contour.setIdentity();
    result.contour -= result.lag;
    return result;
}

/**
 * @brief Convert Log(T_ref^-1 T) into independent Cartesian-position and SO(3)
 *        residuals [p_ref^T(p-p_ref); Log(R_ref^T R)].
 */
inline Eigen::Matrix<double,6,1>
rmpcc_decoupled_error(const Eigen::Matrix<double,6,1> &se3Error)
{
    const Eigen::Matrix4d relative =
        RobotLibrary::Math::se3_exponential(se3Error);
    Eigen::Matrix<double,6,1> result;
    result.head<3>() = relative.block<3,1>(0,3);
    result.tail<3>() = RobotLibrary::Math::so3_logarithm(
        relative.block<3,3>(0,0));
    return result;
}

/** Numerical differential of rmpcc_decoupled_error at the nominal error. */
inline Eigen::Matrix<double,6,6>
rmpcc_decoupled_error_jacobian(
    const Eigen::Matrix<double,6,1> &se3Error,
    const double finiteDifferenceStep)
{
    Eigen::Matrix<double,6,6> jacobian;
    for(int column = 0; column < 6; ++column)
    {
        Eigen::Matrix<double,6,1> plus = se3Error;
        Eigen::Matrix<double,6,1> minus = se3Error;
        plus(column) += finiteDifferenceStep;
        minus(column) -= finiteDifferenceStep;
        jacobian.col(column) =
            (rmpcc_decoupled_error(plus) - rmpcc_decoupled_error(minus))
            / (2.0 * finiteDifferenceStep);
    }
    return jacobian;
}

/**
 * @brief MPCC-equivalent position contour/lag plus independent SO(3) weight.
 *        The angular tangent is intentionally absent from this projector, so
 *        orientation tracking cannot create a hidden progress penalty.
 */
inline Eigen::Matrix<double,6,6>
rmpcc_decoupled_cost_weight(
    const Eigen::Matrix<double,6,1> &referenceBodyTangent,
    const Eigen::Matrix<double,6,6> &contourWeight,
    const Eigen::Matrix<double,6,6> &lagWeight,
    const double regularization)
{
    Eigen::Vector3d direction = referenceBodyTangent.head<3>();
    const double norm = direction.norm();
    if(norm > regularization)
    {
        direction /= norm;
    }
    else
    {
        direction = Eigen::Vector3d::UnitX();
    }
    const Eigen::Matrix3d lag = direction * direction.transpose();
    const Eigen::Matrix3d contour = Eigen::Matrix3d::Identity() - lag;
    Eigen::Matrix<double,6,6> result = Eigen::Matrix<double,6,6>::Zero();
    result.block<3,3>(0,0) =
        contour.transpose() * contourWeight.block<3,3>(0,0) * contour
        + lag.transpose() * lagWeight.block<3,3>(0,0) * lag;
    result.block<3,3>(3,3) = contourWeight.block<3,3>(3,3);
    return result;
}

} } // namespace RobotLibrary::Control

#endif
