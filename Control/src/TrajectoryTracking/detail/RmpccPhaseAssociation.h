/**
 * @file TrajectoryTracking/detail/RmpccPhaseAssociation.h
 * @brief The single phase-association implementation shared by every RMPCC
 *        residual family.
 *
 * Both the local residual family (RmpccPhaseResidual.h) and the associated
 * residual family (RmpccAssociatedPhase.h) need the same quantity: the scalar
 * progress shift delta_s that best explains the current pose error. Keeping two
 * copies of that formula is what allowed the position-only task association to
 * survive in the associated family after it had been repaired in the local one,
 * so it is defined exactly once here and both families call it.
 *
 * Three definitions are always evaluated so that mechanism diagnostics can
 * compare them; the caller selects which one drives the residual.
 *
 *   MetricScrew      delta_s = (g^T M eta) / (g^T M g)
 *   TaskPointXYZ     delta_s = (t_p^T e_p) / (t_p^T t_p)
 *   TaskPoseFeature  delta_s = (h^T dy) / (h^T h),
 *                    y(T) = [p; l_R R e_1; l_R R e_2]
 *
 * TaskPointXYZ divides by the position tangent alone, so on a rotation-rich
 * path it under-normalises by exactly the angular-to-linear tangent ratio.
 * TaskPoseFeature restores the rotational content of the denominator.
 */

#ifndef RMPCC_PHASE_ASSOCIATION_H
#define RMPCC_PHASE_ASSOCIATION_H

#include <Control/TrajectoryTracking/RmpccTypes.h>
#include <Math/MathFunctions.h>

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace RobotLibrary { namespace Control {

struct RmpccPhaseAssociationState
{
    // Shared geometry, exposed because callers scale scalar lags with it.
    Eigen::Vector3d taskTangent = Eigen::Vector3d::Zero();   ///< R_ref * tau_v = dp_ref/ds
    Eigen::Vector3d taskError   = Eigen::Vector3d::Zero();   ///< p_actual - p_ref

    // All three definitions, always evaluated for mechanism diagnostics.
    double metricCorrection = 0.0;
    double metricDenominator = 0.0;
    double taskCorrection = 0.0;
    double taskDenominator = 0.0;
    bool   taskObservable = false;
    double featureCorrection = 0.0;
    double featureDenominator = 0.0;
    bool   featureObservable = false;

    // The definition selected by phaseAssociation.
    double correction = 0.0;
    double denominator = 0.0;
    bool   observable = false;
};

/**
 * @param error          eta = Log(D(s)^-1 X)
 * @param errorTangent   g = J_r(eta)^-1 Ad(Exp(eta)^-1) tau(s)
 * @param reference      D(s)
 * @param actual         X = D(s) Exp(eta)
 * @param referenceTangent tau(s), the reference body tangent
 * @param rotationLength l_R, required positive only for TaskPoseFeature
 */
inline RmpccPhaseAssociationState
rmpcc_evaluate_phase_association(
    const Eigen::Vector<double,6> &error,
    const Eigen::Vector<double,6> &errorTangent,
    const Eigen::Matrix4d &reference,
    const Eigen::Matrix4d &actual,
    const Eigen::Vector<double,6> &referenceTangent,
    const Eigen::Matrix<double,6,6> &metric,
    const RmpccPhaseAssociation phaseAssociation,
    const double metricRegularization,
    const double phaseDenominatorTolerance,
    const double rotationLength)
{
    RmpccPhaseAssociationState result;

    result.metricDenominator =
        (errorTangent.transpose() * metric * errorTangent)(0);
    result.metricCorrection =
        (errorTangent.transpose() * metric * error)(0)
        / std::max(result.metricDenominator, metricRegularization);

    result.taskError = actual.block<3,1>(0,3) - reference.block<3,1>(0,3);
    result.taskTangent =
        reference.block<3,3>(0,0) * referenceTangent.head<3>();
    result.taskDenominator = result.taskTangent.squaredNorm();
    result.taskObservable = std::isfinite(result.taskDenominator)
        && result.taskDenominator > phaseDenominatorTolerance;
    if(result.taskObservable)
    {
        result.taskCorrection =
            result.taskTangent.dot(result.taskError) / result.taskDenominator;
    }

    if(phaseAssociation == RmpccPhaseAssociation::TaskPoseFeature
       && (not std::isfinite(rotationLength) or rotationLength <= 0.0))
    {
        throw std::invalid_argument(
            "[ERROR] [RMPCC PHASE RESIDUAL] TaskPoseFeature rotation length must be positive.");
    }
    const Eigen::Matrix3d referenceRotation = reference.block<3,3>(0,0);
    const Eigen::Matrix3d actualRotation = actual.block<3,3>(0,0);
    const Eigen::Vector3d angularTangent = referenceTangent.tail<3>();
    Eigen::Vector<double,9> featureTangent = Eigen::Vector<double,9>::Zero();
    Eigen::Vector<double,9> featureError = Eigen::Vector<double,9>::Zero();
    featureTangent.head<3>() = result.taskTangent;
    featureError.head<3>() = result.taskError;
    for(int axis = 0; axis < 2; ++axis)
    {
        const Eigen::Vector3d basis = Eigen::Vector3d::Unit(axis);
        featureTangent.segment<3>(3 + 3 * axis) =
            rotationLength * referenceRotation * angularTangent.cross(basis);
        featureError.segment<3>(3 + 3 * axis) =
            rotationLength
            * (actualRotation.col(axis) - referenceRotation.col(axis));
    }
    result.featureDenominator = featureTangent.squaredNorm();
    result.featureObservable = std::isfinite(result.featureDenominator)
        && result.featureDenominator > phaseDenominatorTolerance;
    if(result.featureObservable)
    {
        result.featureCorrection =
            featureTangent.dot(featureError) / result.featureDenominator;
    }

    switch(phaseAssociation)
    {
        case RmpccPhaseAssociation::MetricScrew:
            result.correction = result.metricCorrection;
            result.denominator = result.metricDenominator;
            result.observable = std::isfinite(result.metricCorrection);
            break;
        case RmpccPhaseAssociation::TaskPointXYZ:
            result.correction = result.taskCorrection;
            result.denominator = result.taskDenominator;
            result.observable = result.taskObservable;
            break;
        case RmpccPhaseAssociation::TaskPoseFeature:
            result.correction = result.featureCorrection;
            result.denominator = result.featureDenominator;
            result.observable = result.featureObservable;
            break;
    }
    return result;
}

} } // namespace RobotLibrary::Control

#endif
