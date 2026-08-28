/**
 * @file TrajectoryTracking/detail/RmpccPhaseResidual.h
 * @brief Pure task-phase and lag-representation residual geometry for RMPCC.
 */

#ifndef RMPCC_PHASE_RESIDUAL_H
#define RMPCC_PHASE_RESIDUAL_H

#include "RmpccPhaseAssociation.h"
#include "RmpccPrediction.h"

#include <Control/TrajectoryTracking/RmpccTypes.h>
#include <Math/MathFunctions.h>

#include <Eigen/Core>

#include <algorithm>
#include <array>
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

inline RmpccPhaseResiduals
rmpcc_metric_phase_residuals_from_tangent(
    const Eigen::Vector<double,6> &eta,
    const Eigen::Vector<double,6> &tau,
    const Eigen::Matrix<double,6,6> &metric,
    const double metricRegularization)
{
    const Eigen::Vector<double,6> g =
        rmpcc_error_coordinate_path_tangent(eta, tau);
    RmpccPhaseResiduals result;
    result.metricPhaseDenominator = (g.transpose() * metric * g)(0);
    result.metricPhaseCorrection =
        (g.transpose() * metric * eta)(0)
        / std::max(result.metricPhaseDenominator, metricRegularization);
    result.phaseCorrection = result.metricPhaseCorrection;
    result.phaseDenominator = result.metricPhaseDenominator;
    result.phaseObservable = std::isfinite(result.phaseCorrection);
    if(!result.phaseObservable)
    {
        result.contour = eta;
        return result;
    }
    result.vectorLag = g * result.phaseCorrection;
    result.scalarLag = tau.head<3>().norm() * result.phaseCorrection;
    result.contour = eta - result.vectorLag;
    return result;
}

inline void
rmpcc_linearize_metric_phase_residuals_eta(
    const Eigen::Vector<double,6> &eta,
    const Eigen::Vector<double,6> &tau,
    const Eigen::Matrix<double,6,6> &metric,
    const double metricRegularization,
    RmpccPhaseResidualLinearization &result)
{
    using Matrix6d = Eigen::Matrix<double,6,6>;
    using Vector6d = Eigen::Vector<double,6>;
    const Matrix6d ad = RobotLibrary::Math::se3_adjoint_matrix(eta);
    constexpr std::array<double,13> coefficients = {
        1.0, -0.5, 1.0 / 12.0, 0.0, -1.0 / 720.0, 0.0,
        1.0 / 30240.0, 0.0, -1.0 / 1209600.0, 0.0,
        1.0 / 47900160.0, 0.0, -691.0 / 1307674368000.0};

    std::array<Vector6d,13> powers;
    powers[0] = tau;
    Vector6d g = coefficients[0] * powers[0];
    for(int order = 1; order <= 12; ++order)
    {
        powers[static_cast<size_t>(order)] =
            ad * powers[static_cast<size_t>(order - 1)];
        if(coefficients[static_cast<size_t>(order)] != 0.0)
        {
            g += coefficients[static_cast<size_t>(order)]
                * powers[static_cast<size_t>(order)];
        }
    }

    Matrix6d dg;
    for(int column = 0; column < 6; ++column)
    {
        Vector6d basis = Vector6d::Zero();
        basis(column) = 1.0;
        const Matrix6d dad = RobotLibrary::Math::se3_adjoint_matrix(basis);
        Vector6d derivative = Vector6d::Zero();
        Vector6d dgColumn = Vector6d::Zero();
        for(int order = 1; order <= 12; ++order)
        {
            derivative = dad * powers[static_cast<size_t>(order - 1)]
                + ad * derivative;
            if(coefficients[static_cast<size_t>(order)] != 0.0)
            {
                dgColumn += coefficients[static_cast<size_t>(order)] * derivative;
            }
        }
        dg.col(column) = dgColumn;
    }

    const double rawDenominator = (g.transpose() * metric * g)(0);
    const double denominator = std::max(rawDenominator, metricRegularization);
    const double numerator = (g.transpose() * metric * eta)(0);
    const double alpha = numerator / denominator;
    for(int column = 0; column < 6; ++column)
    {
        Vector6d basis = Vector6d::Zero();
        basis(column) = 1.0;
        const double dNumerator =
            (dg.col(column).transpose() * metric * eta)(0)
            + (g.transpose() * metric * basis)(0);
        const double dDenominator = rawDenominator > metricRegularization
            ? 2.0 * (dg.col(column).transpose() * metric * g)(0) : 0.0;
        const double dAlpha =
            (dNumerator * denominator - numerator * dDenominator)
            / (denominator * denominator);
        const Vector6d dLag = dg.col(column) * alpha + g * dAlpha;
        result.vectorLagJacobian.col(column) = dLag;
        result.contourJacobian.col(column) = basis - dLag;
        result.scalarLagJacobian(0, column) = tau.head<3>().norm() * dAlpha;
    }
}

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
    const Eigen::Vector<double,6> tau = referenceTangent(state(6));
    const Eigen::Vector<double,6> g =
        rmpcc_error_coordinate_path_tangent(eta, tau);

    // The final MetricScrew residual depends only on eta and the transported
    // body tangent.  Constructing reference/actual poses here was solely for
    // the TaskPointXYZ/TaskPoseFeature research alternatives and multiplied
    // that diagnostic work across every finite-difference column.
    if(phaseAssociation == RmpccPhaseAssociation::MetricScrew)
    {
        return rmpcc_metric_phase_residuals_from_tangent(
            eta, tau, metric, metricRegularization);
    }

    const Eigen::Matrix4d reference = referenceTransform(state(6));
    const Eigen::Matrix4d actual =
        reference * RobotLibrary::Math::se3_exponential(eta);

    const RmpccPhaseAssociationState phase = rmpcc_evaluate_phase_association(
        eta, g, reference, actual, tau, metric, phaseAssociation,
        metricRegularization, phaseDenominatorTolerance,
        taskPoseFeatureRotationLength);

    RmpccPhaseResiduals result;
    result.metricPhaseDenominator = phase.metricDenominator;
    result.metricPhaseCorrection = phase.metricCorrection;
    result.taskPhaseDenominator = phase.taskDenominator;
    result.taskPhaseObservable = phase.taskObservable;
    result.taskPhaseCorrection = phase.taskCorrection;
    result.taskPoseFeaturePhaseDenominator = phase.featureDenominator;
    result.taskPoseFeaturePhaseObservable = phase.featureObservable;
    result.taskPoseFeaturePhaseCorrection = phase.featureCorrection;

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

    result.phaseCorrection = phase.correction;
    result.phaseDenominator = phase.denominator;
    result.phaseObservable = phase.observable;
    // ScalarTaskDistance semantics are deliberately unchanged: the scalar lag
    // is the physical task-point along-path distance |t_p| delta_s, in metres,
    // whichever association produced delta_s.
    result.vectorLag = g * result.phaseCorrection;
    result.scalarLag = phase.taskTangent.norm() * result.phaseCorrection;
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
    if(phaseAssociation == RmpccPhaseAssociation::MetricScrew)
    {
        const Eigen::Vector<double,6> nominalTangent = tangent(state(6));
        if(state.head<6>().tail<3>().norm() <= 0.5)
        {
            rmpcc_linearize_metric_phase_residuals_eta(
                state.head<6>(), nominalTangent, metric,
                metricRegularization, result);
        }
        else
        {
            for(int column = 0; column < 6; ++column)
            {
                Eigen::Vector<double,6> plus = state.head<6>();
                Eigen::Vector<double,6> minus = state.head<6>();
                plus(column) += finiteDifferenceStep;
                minus(column) -= finiteDifferenceStep;
                const RmpccPhaseResiduals plusResidual =
                    rmpcc_metric_phase_residuals_from_tangent(
                        plus, nominalTangent, metric, metricRegularization);
                const RmpccPhaseResiduals minusResidual =
                    rmpcc_metric_phase_residuals_from_tangent(
                        minus, nominalTangent, metric, metricRegularization);
                const double etaDenominator = 2.0 * finiteDifferenceStep;
                result.contourJacobian.col(column) =
                    (plusResidual.contour - minusResidual.contour) / etaDenominator;
                result.vectorLagJacobian.col(column) =
                    (plusResidual.vectorLag - minusResidual.vectorLag) / etaDenominator;
                result.scalarLagJacobian(0, column) =
                    (plusResidual.scalarLag - minusResidual.scalarLag)
                    / etaDenominator;
            }
        }
        RmpccStateVector plus = state;
        RmpccStateVector minus = state;
        plus(6) = std::clamp(plus(6) + finiteDifferenceStep, 0.0, 1.0);
        minus(6) = std::clamp(minus(6) - finiteDifferenceStep, 0.0, 1.0);
        const double denominator = plus(6) - minus(6);
        if(std::abs(denominator) > 1e-15)
        {
            const RmpccPhaseResiduals plusResidual = rmpcc_phase_residuals(
                plus, metric, phaseAssociation, metricRegularization,
                phaseDenominatorTolerance, taskPoseFeatureRotationLength,
                transform, tangent);
            const RmpccPhaseResiduals minusResidual = rmpcc_phase_residuals(
                minus, metric, phaseAssociation, metricRegularization,
                phaseDenominatorTolerance, taskPoseFeatureRotationLength,
                transform, tangent);
            result.contourJacobian.col(6) =
                (plusResidual.contour - minusResidual.contour) / denominator;
            result.vectorLagJacobian.col(6) =
                (plusResidual.vectorLag - minusResidual.vectorLag) / denominator;
            result.scalarLagJacobian(0, 6) =
                (plusResidual.scalarLag - minusResidual.scalarLag) / denominator;
        }
        return result;
    }
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
