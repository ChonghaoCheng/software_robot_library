/**
 * @file TrajectoryTracking/detail/RmpccAssociatedPhase.h
 * @brief Shared TaskPoint-associated pose context and pose-path arc geometry.
 */

#ifndef RMPCC_ASSOCIATED_PHASE_H
#define RMPCC_ASSOCIATED_PHASE_H

#include "RmpccCostGeometry.h"
#include "RmpccPhaseAssociation.h"
#include "RmpccPrediction.h"

#include <Control/TrajectoryTracking/RmpccTypes.h>
#include <Math/MathFunctions.h>

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace RobotLibrary { namespace Control {

class RmpccPoseArcTable
{
    public:
        template<typename ReferenceTangentFunction>
        void build(const Eigen::Matrix<double,6,6> &lagWeight,
                   ReferenceTangentFunction &&referenceTangent,
                   const std::size_t sampleCount = 4097)
        {
            if(sampleCount < 2)
            {
                throw std::invalid_argument(
                    "[ERROR] [RMPCC POSE ARC] At least two samples are required.");
            }
            const Eigen::Vector3d translationWeights =
                lagWeight.diagonal().head<3>();
            const double tolerance =
                1e-12 * std::max(1.0, std::abs(translationWeights(0)));
            if(not std::isfinite(translationWeights(0))
               or translationWeights(0) <= 0.0
               or (translationWeights.array() - translationWeights(0))
                      .abs().maxCoeff() > tolerance)
            {
                throw std::invalid_argument(
                    "[ERROR] [RMPCC POSE ARC] Equal positive translational lag weights are required.");
            }

            _qLag = translationWeights(0);
            _metric = lagWeight / _qLag;
            _sampleCount = sampleCount;
            _arc.assign(sampleCount, 0.0);
            _density.assign(sampleCount, 0.0);
            _minimumDensity = std::numeric_limits<double>::infinity();
            _maximumDensity = 0.0;
            const double step = 1.0 / static_cast<double>(sampleCount - 1);
            for(std::size_t index = 0; index < sampleCount; ++index)
            {
                const double progress = step * static_cast<double>(index);
                const Eigen::Vector<double,6> tangent =
                    referenceTangent(progress);
                const double quadratic = (tangent.transpose() * _metric * tangent)(0);
                _density[index] = std::sqrt(std::max(0.0, quadratic));
                _minimumDensity = std::min(_minimumDensity, _density[index]);
                _maximumDensity = std::max(_maximumDensity, _density[index]);
                if(index > 0)
                {
                    _arc[index] = _arc[index - 1]
                        + 0.5 * step * (_density[index - 1] + _density[index]);
                }
            }
        }

        double value(const double progress) const
        {
            if(_arc.empty())
            {
                throw std::logic_error("[ERROR] [RMPCC POSE ARC] Table is not built.");
            }
            const double clamped = std::clamp(progress, 0.0, 1.0);
            const double coordinate =
                clamped * static_cast<double>(_sampleCount - 1);
            const std::size_t lower = static_cast<std::size_t>(coordinate);
            if(lower >= _sampleCount - 1)
            {
                return _arc.back();
            }
            const double fraction = coordinate - static_cast<double>(lower);
            return (1.0 - fraction) * _arc[lower] + fraction * _arc[lower + 1];
        }

        double lag(const double internalProgress,
                   const double associatedProgress) const
        {
            return value(associatedProgress) - value(internalProgress);
        }

        std::size_t sample_count() const { return _sampleCount; }
        double q_lag() const { return _qLag; }
        double minimum_density() const { return _minimumDensity; }
        double maximum_density() const { return _maximumDensity; }
        double total_arc_length() const { return _arc.empty() ? 0.0 : _arc.back(); }
        const Eigen::Matrix<double,6,6>& metric() const { return _metric; }

    private:
        std::size_t _sampleCount = 0;
        double _qLag = 0.0;
        double _minimumDensity = 0.0;
        double _maximumDensity = 0.0;
        Eigen::Matrix<double,6,6> _metric =
            Eigen::Matrix<double,6,6>::Identity();
        std::vector<double> _arc;
        std::vector<double> _density;
};

struct RmpccAssociatedPhaseContext
{
    double internalProgress = 0.0;
    double phaseCorrection = 0.0;
    double associatedProgress = 0.0;
    Eigen::Matrix4d referenceAtInternalProgress = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d referenceAtAssociatedProgress = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d predictedActualPose = Eigen::Matrix4d::Identity();
    Eigen::Vector<double,6> tangentAtInternalProgress =
        Eigen::Vector<double,6>::Zero();
    Eigen::Vector<double,6> tangentAtAssociatedProgress =
        Eigen::Vector<double,6>::Zero();
    Eigen::Vector<double,6> errorTangentAtInternalProgress =
        Eigen::Vector<double,6>::Zero();
    double phaseDenominator = 0.0;
    double positionLag = 0.0;
    double poseArcLag = 0.0;
    bool observable = false;
};

template<typename ReferenceTransformFunction, typename ReferenceTangentFunction>
RmpccAssociatedPhaseContext
rmpcc_associated_phase_context(
    const RmpccStateVector &state,
    const Eigen::Matrix<double,6,6> &metric,
    const RmpccPhaseAssociation phaseAssociation,
    const double metricRegularization,
    const double phaseDenominatorTolerance,
    const double taskPoseFeatureRotationLength,
    const RmpccPoseArcTable &poseArc,
    ReferenceTransformFunction &&referenceTransform,
    ReferenceTangentFunction &&referenceTangent)
{
    RmpccAssociatedPhaseContext result;
    result.internalProgress = state(6);
    result.referenceAtInternalProgress = referenceTransform(state(6));
    result.tangentAtInternalProgress = referenceTangent(state(6));
    result.predictedActualPose = result.referenceAtInternalProgress
        * RobotLibrary::Math::se3_exponential(state.head<6>());

    result.errorTangentAtInternalProgress =
        rmpcc_error_coordinate_path_tangent(
            state.head<6>(), result.tangentAtInternalProgress);

    const RmpccPhaseAssociationState phase = rmpcc_evaluate_phase_association(
        state.head<6>(), result.errorTangentAtInternalProgress,
        result.referenceAtInternalProgress, result.predictedActualPose,
        result.tangentAtInternalProgress, metric, phaseAssociation,
        metricRegularization, phaseDenominatorTolerance,
        taskPoseFeatureRotationLength);
    const Eigen::Vector3d taskTangent = phase.taskTangent;
    result.phaseDenominator = phase.denominator;
    result.observable = phase.observable;
    result.phaseCorrection = phase.correction;
    result.associatedProgress = std::clamp(
        result.internalProgress + result.phaseCorrection, 0.0, 1.0);
    result.referenceAtAssociatedProgress =
        referenceTransform(result.associatedProgress);
    result.tangentAtAssociatedProgress =
        referenceTangent(result.associatedProgress);
    result.positionLag = taskTangent.norm() * result.phaseCorrection;
    result.poseArcLag = poseArc.lag(
        result.internalProgress, result.associatedProgress);
    return result;
}

struct RmpccAssociatedResiduals
{
    RmpccAssociatedPhaseContext context;
    Eigen::Vector<double,6> contour = Eigen::Vector<double,6>::Zero();
    Eigen::Vector<double,6> vectorLag = Eigen::Vector<double,6>::Zero();
    double scalarPoseArcLag = 0.0;
};

struct RmpccAssociatedResidualLinearization
{
    RmpccAssociatedResiduals residual;
    Eigen::Matrix<double,6,7> contourJacobian =
        Eigen::Matrix<double,6,7>::Zero();
    Eigen::Matrix<double,6,7> vectorLagJacobian =
        Eigen::Matrix<double,6,7>::Zero();
    Eigen::Matrix<double,1,7> scalarPoseArcLagJacobian =
        Eigen::Matrix<double,1,7>::Zero();
};

template<typename ReferenceTransformFunction, typename ReferenceTangentFunction>
RmpccAssociatedResiduals
rmpcc_associated_residuals(
    const RmpccStateVector &state,
    const RmpccContourResidualGeometry contourGeometry,
    const Eigen::Matrix<double,6,6> &metric,
    const RmpccPhaseAssociation phaseAssociation,
    const double metricRegularization,
    const double phaseDenominatorTolerance,
    const double taskPoseFeatureRotationLength,
    const RmpccPoseArcTable &poseArc,
    ReferenceTransformFunction &&referenceTransform,
    ReferenceTangentFunction &&referenceTangent)
{
    RmpccAssociatedResiduals result;
    result.context = rmpcc_associated_phase_context(
        state, metric, phaseAssociation, metricRegularization,
        phaseDenominatorTolerance, taskPoseFeatureRotationLength, poseArc,
        referenceTransform, referenceTangent);
    const Eigen::Matrix4d relative =
        RobotLibrary::Math::se3_inverse(
            result.context.referenceAtAssociatedProgress)
        * result.context.predictedActualPose;
    const Eigen::Vector<double,6> exactGroup =
        RobotLibrary::Math::se3_logarithm(relative);
    result.contour = contourGeometry
            == RmpccContourResidualGeometry::AssociatedDecoupledCartesianSO3
        ? rmpcc_decoupled_error(exactGroup) : exactGroup;
    result.vectorLag = result.context.errorTangentAtInternalProgress
        * result.context.phaseCorrection;
    result.scalarPoseArcLag = result.context.poseArcLag;
    return result;
}

template<typename ReferenceTransformFunction, typename ReferenceTangentFunction>
RmpccAssociatedResidualLinearization
rmpcc_linearize_associated_residuals(
    const RmpccStateVector &state,
    const RmpccContourResidualGeometry contourGeometry,
    const Eigen::Matrix<double,6,6> &metric,
    const RmpccPhaseAssociation phaseAssociation,
    const double metricRegularization,
    const double phaseDenominatorTolerance,
    const double taskPoseFeatureRotationLength,
    const double finiteDifferenceStep,
    const RmpccPoseArcTable &poseArc,
    ReferenceTransformFunction &&referenceTransform,
    ReferenceTangentFunction &&referenceTangent)
{
    if(not std::isfinite(finiteDifferenceStep) or finiteDifferenceStep <= 0.0)
    {
        throw std::invalid_argument(
            "[ERROR] [RMPCC ASSOCIATED PHASE] finiteDifferenceStep must be positive.");
    }
    auto &&transform = referenceTransform;
    auto &&tangent = referenceTangent;
    RmpccAssociatedResidualLinearization result;
    result.residual = rmpcc_associated_residuals(
        state, contourGeometry, metric, phaseAssociation,
        metricRegularization, phaseDenominatorTolerance,
        taskPoseFeatureRotationLength, poseArc, transform, tangent);
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
        const RmpccAssociatedResiduals plusResidual =
            rmpcc_associated_residuals(
                plus, contourGeometry, metric, phaseAssociation,
                metricRegularization, phaseDenominatorTolerance,
                taskPoseFeatureRotationLength, poseArc, transform, tangent);
        const RmpccAssociatedResiduals minusResidual =
            rmpcc_associated_residuals(
                minus, contourGeometry, metric, phaseAssociation,
                metricRegularization, phaseDenominatorTolerance,
                taskPoseFeatureRotationLength, poseArc, transform, tangent);
        result.contourJacobian.col(column) =
            (plusResidual.contour - minusResidual.contour) / denominator;
        result.vectorLagJacobian.col(column) =
            (plusResidual.vectorLag - minusResidual.vectorLag) / denominator;
        result.scalarPoseArcLagJacobian(0, column) =
            (plusResidual.scalarPoseArcLag - minusResidual.scalarPoseArcLag)
            / denominator;
    }
    return result;
}

/** Legacy TaskPointXYZ-only overloads; behaviour is bit-identical to before. */
template<typename ReferenceTransformFunction, typename ReferenceTangentFunction>
RmpccAssociatedResiduals
rmpcc_associated_residuals(
    const RmpccStateVector &state,
    const RmpccContourResidualGeometry contourGeometry,
    const double phaseDenominatorTolerance,
    const RmpccPoseArcTable &poseArc,
    ReferenceTransformFunction &&referenceTransform,
    ReferenceTangentFunction &&referenceTangent)
{
    return rmpcc_associated_residuals(
        state, contourGeometry, Eigen::Matrix<double,6,6>::Identity(),
        RmpccPhaseAssociation::TaskPointXYZ, 0.0, phaseDenominatorTolerance,
        0.0, poseArc,
        std::forward<ReferenceTransformFunction>(referenceTransform),
        std::forward<ReferenceTangentFunction>(referenceTangent));
}

template<typename ReferenceTransformFunction, typename ReferenceTangentFunction>
RmpccAssociatedResidualLinearization
rmpcc_linearize_associated_residuals(
    const RmpccStateVector &state,
    const RmpccContourResidualGeometry contourGeometry,
    const double phaseDenominatorTolerance,
    const double finiteDifferenceStep,
    const RmpccPoseArcTable &poseArc,
    ReferenceTransformFunction &&referenceTransform,
    ReferenceTangentFunction &&referenceTangent)
{
    return rmpcc_linearize_associated_residuals(
        state, contourGeometry, Eigen::Matrix<double,6,6>::Identity(),
        RmpccPhaseAssociation::TaskPointXYZ, 0.0, phaseDenominatorTolerance,
        0.0, finiteDifferenceStep, poseArc,
        std::forward<ReferenceTransformFunction>(referenceTransform),
        std::forward<ReferenceTangentFunction>(referenceTangent));
}

} } // namespace RobotLibrary::Control

#endif
