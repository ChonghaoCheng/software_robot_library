/**
 * @file TrajectoryTracking/detail/RmpccPhaseResidual.h
 * @brief Pure task-phase and lag-representation residual geometry for RMPCC.
 */

#ifndef RMPCC_PHASE_RESIDUAL_H
#define RMPCC_PHASE_RESIDUAL_H

#include "RmpccPrediction.h"

#include <Control/TrajectoryTracking/RmpccTypes.h>
#include <Math/MathFunctions.h>

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace RobotLibrary { namespace Control {

struct RmpccPhaseResiduals
{
    Eigen::Vector<double,6> contour = Eigen::Vector<double,6>::Zero();
    Eigen::Vector<double,6> vectorLag = Eigen::Vector<double,6>::Zero();
    double scalarLag = 0.0;
    double phaseCorrection = 0.0;
    double phaseDenominator = 0.0;
    bool phaseObservable = false;

    // Both phase definitions are always evaluated for mechanism diagnostics.
    double metricPhaseCorrection = 0.0;
    double metricPhaseDenominator = 0.0;
    double taskPhaseCorrection = 0.0;
    double taskPhaseDenominator = 0.0;
    bool taskPhaseObservable = false;
    double taskPoseFeaturePhaseCorrection = 0.0;
    double taskPoseFeaturePhaseDenominator = 0.0;
    bool taskPoseFeaturePhaseObservable = false;
};

struct RmpccPhaseResidualLinearization
{
    RmpccPhaseResiduals residual;
    Eigen::Matrix<double,6,7> contourJacobian =
        Eigen::Matrix<double,6,7>::Zero();
    Eigen::Matrix<double,6,7> vectorLagJacobian =
        Eigen::Matrix<double,6,7>::Zero();
    Eigen::Matrix<double,1,7> scalarLagJacobian =
        Eigen::Matrix<double,1,7>::Zero();
};

template<typename ReferenceTransformFunction, typename ReferenceTangentFunction>
RmpccPhaseResiduals
rmpcc_phase_residuals(
    const RmpccStateVector &state,
    const Eigen::Matrix<double,6,6> &metric,
    const RmpccPhaseAssociation phaseAssociation,
    const double metricRegularization,
    const double phaseDenominatorTolerance,
    const double taskPoseFeatureRotationLength,
    ReferenceTransformFunction &&referenceTransform,
    ReferenceTangentFunction &&referenceTangent)
{
    const Eigen::Vector<double,6> eta = state.head<6>();
    const Eigen::Matrix4d reference = referenceTransform(state(6));
    const Eigen::Vector<double,6> tau = referenceTangent(state(6));
    const Eigen::Vector<double,6> g =
        rmpcc_error_coordinate_path_tangent(eta, tau);

    RmpccPhaseResiduals result;
    result.metricPhaseDenominator = (g.transpose() * metric * g)(0);
    const double metricDenominator =
        std::max(result.metricPhaseDenominator, metricRegularization);
    result.metricPhaseCorrection =
        (g.transpose() * metric * eta)(0) / metricDenominator;

    const Eigen::Matrix4d actual =
        reference * RobotLibrary::Math::se3_exponential(eta);
    const Eigen::Vector3d taskError =
        actual.block<3,1>(0,3) - reference.block<3,1>(0,3);
    const Eigen::Vector3d taskTangent =
        reference.block<3,3>(0,0) * tau.head<3>();
    result.taskPhaseDenominator = taskTangent.squaredNorm();
    result.taskPhaseObservable =
        std::isfinite(result.taskPhaseDenominator)
        && result.taskPhaseDenominator > phaseDenominatorTolerance;
    if(result.taskPhaseObservable)
    {
        result.taskPhaseCorrection =
            taskTangent.dot(taskError) / result.taskPhaseDenominator;
    }

    if(phaseAssociation == RmpccPhaseAssociation::TaskPoseFeature
       && (not std::isfinite(taskPoseFeatureRotationLength)
           or taskPoseFeatureRotationLength <= 0.0))
    {
        throw std::invalid_argument(
            "[ERROR] [RMPCC PHASE RESIDUAL] TaskPoseFeature rotation length must be positive.");
    }
    const Eigen::Matrix3d referenceRotation = reference.block<3,3>(0,0);
    const Eigen::Matrix3d actualRotation = actual.block<3,3>(0,0);
    const Eigen::Vector3d angularTangent = tau.tail<3>();
    Eigen::Vector<double,9> featureTangent =
        Eigen::Vector<double,9>::Zero();
    Eigen::Vector<double,9> featureError =
        Eigen::Vector<double,9>::Zero();
    featureTangent.head<3>() = taskTangent;
    featureError.head<3>() = taskError;
    for(int axis = 0; axis < 2; ++axis)
    {
        const Eigen::Vector3d basis = Eigen::Vector3d::Unit(axis);
        featureTangent.segment<3>(3 + 3 * axis) =
            taskPoseFeatureRotationLength * referenceRotation
            * angularTangent.cross(basis);
        featureError.segment<3>(3 + 3 * axis) =
            taskPoseFeatureRotationLength
            * (actualRotation.col(axis) - referenceRotation.col(axis));
    }
    result.taskPoseFeaturePhaseDenominator = featureTangent.squaredNorm();
    result.taskPoseFeaturePhaseObservable =
        std::isfinite(result.taskPoseFeaturePhaseDenominator)
        && result.taskPoseFeaturePhaseDenominator > phaseDenominatorTolerance;
    if(result.taskPoseFeaturePhaseObservable)
    {
        result.taskPoseFeaturePhaseCorrection =
            featureTangent.dot(featureError)
            / result.taskPoseFeaturePhaseDenominator;
    }

    const bool selectedObservable =
        phaseAssociation == RmpccPhaseAssociation::MetricScrew
        or (phaseAssociation == RmpccPhaseAssociation::TaskPointXYZ
            && result.taskPhaseObservable)
        or (phaseAssociation == RmpccPhaseAssociation::TaskPoseFeature
            && result.taskPoseFeaturePhaseObservable);
    if(not selectedObservable)
    {
        result.contour = eta;
        result.phaseDenominator =
            phaseAssociation == RmpccPhaseAssociation::TaskPointXYZ
            ? result.taskPhaseDenominator
            : result.taskPoseFeaturePhaseDenominator;
        return result;
    }

    if(phaseAssociation == RmpccPhaseAssociation::MetricScrew)
    {
        result.phaseCorrection = result.metricPhaseCorrection;
        result.phaseDenominator = result.metricPhaseDenominator;
        result.phaseObservable = std::isfinite(result.metricPhaseCorrection);
    }
    else if(phaseAssociation == RmpccPhaseAssociation::TaskPointXYZ)
    {
        result.phaseCorrection = result.taskPhaseCorrection;
        result.phaseDenominator = result.taskPhaseDenominator;
        result.phaseObservable = result.taskPhaseObservable;
    }
    else
    {
        result.phaseCorrection = result.taskPoseFeaturePhaseCorrection;
        result.phaseDenominator = result.taskPoseFeaturePhaseDenominator;
        result.phaseObservable = result.taskPoseFeaturePhaseObservable;
    }
    result.vectorLag = g * result.phaseCorrection;
    result.scalarLag = taskTangent.norm() * result.phaseCorrection;
    result.contour = eta - result.vectorLag;
    return result;
}

template<typename ReferenceTransformFunction, typename ReferenceTangentFunction>
RmpccPhaseResiduals
rmpcc_phase_residuals(
    const RmpccStateVector &state,
    const Eigen::Matrix<double,6,6> &metric,
    const RmpccPhaseAssociation phaseAssociation,
    const double metricRegularization,
    const double phaseDenominatorTolerance,
    ReferenceTransformFunction &&referenceTransform,
    ReferenceTangentFunction &&referenceTangent)
{
    return rmpcc_phase_residuals(
        state, metric, phaseAssociation, metricRegularization,
        phaseDenominatorTolerance, 0.0,
        std::forward<ReferenceTransformFunction>(referenceTransform),
        std::forward<ReferenceTangentFunction>(referenceTangent));
}

template<typename ReferenceTransformFunction, typename ReferenceTangentFunction>
RmpccPhaseResidualLinearization
rmpcc_linearize_phase_residuals(
    const RmpccStateVector &state,
    const Eigen::Matrix<double,6,6> &metric,
    const RmpccPhaseAssociation phaseAssociation,
    const double metricRegularization,
    const double phaseDenominatorTolerance,
    const double taskPoseFeatureRotationLength,
    const double finiteDifferenceStep,
    ReferenceTransformFunction &&referenceTransform,
    ReferenceTangentFunction &&referenceTangent)
{
    if(not std::isfinite(finiteDifferenceStep) or finiteDifferenceStep <= 0.0)
    {
        throw std::invalid_argument(
            "[ERROR] [RMPCC PHASE RESIDUAL] finiteDifferenceStep must be positive.");
    }

    auto &&transform = referenceTransform;
    auto &&tangent = referenceTangent;
    RmpccPhaseResidualLinearization result;
    result.residual = rmpcc_phase_residuals(
        state, metric, phaseAssociation, metricRegularization,
        phaseDenominatorTolerance, taskPoseFeatureRotationLength,
        transform, tangent);
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
        const RmpccPhaseResiduals plusResidual = rmpcc_phase_residuals(
            plus, metric, phaseAssociation, metricRegularization,
            phaseDenominatorTolerance, taskPoseFeatureRotationLength,
            transform, tangent);
        const RmpccPhaseResiduals minusResidual = rmpcc_phase_residuals(
            minus, metric, phaseAssociation, metricRegularization,
            phaseDenominatorTolerance, taskPoseFeatureRotationLength,
            transform, tangent);
        result.contourJacobian.col(column) =
            (plusResidual.contour - minusResidual.contour) / denominator;
        result.vectorLagJacobian.col(column) =
            (plusResidual.vectorLag - minusResidual.vectorLag) / denominator;
        result.scalarLagJacobian(0, column) =
            (plusResidual.scalarLag - minusResidual.scalarLag) / denominator;
    }
    return result;
}

template<typename ReferenceTransformFunction, typename ReferenceTangentFunction>
RmpccPhaseResidualLinearization
rmpcc_linearize_phase_residuals(
    const RmpccStateVector &state,
    const Eigen::Matrix<double,6,6> &metric,
    const RmpccPhaseAssociation phaseAssociation,
    const double metricRegularization,
    const double phaseDenominatorTolerance,
    const double finiteDifferenceStep,
    ReferenceTransformFunction &&referenceTransform,
    ReferenceTangentFunction &&referenceTangent)
{
    return rmpcc_linearize_phase_residuals(
        state, metric, phaseAssociation, metricRegularization,
        phaseDenominatorTolerance, 0.0, finiteDifferenceStep,
        std::forward<ReferenceTransformFunction>(referenceTransform),
        std::forward<ReferenceTangentFunction>(referenceTangent));
}

} } // namespace RobotLibrary::Control

#endif
