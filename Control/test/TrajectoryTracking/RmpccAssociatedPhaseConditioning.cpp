/**
 * @file TrajectoryTracking/test/RmpccAssociatedPhaseConditioning.cpp
 * @brief F2-02B gate: the associated residual family shares the local family's
 *        phase-association implementation, and is conditioned identically.
 *
 * The associated family (TG/TPA/DT) used to carry its own copy of
 * delta_s = t_p^T e_p / t_p^T t_p, which is why the F2-01 repair of the local
 * family left it untouched. This gate covers the associated family with the
 * same conditioning checks the local family already has, so the two cannot
 * diverge again:
 *
 *   1. local and associated delta_s agree exactly, given the same state;
 *   2. TaskPointXYZ position-to-rotation gain equals the tangent ratio r
 *      (the known-bad baseline, locked in as a characterisation);
 *   3. TaskPoseFeature gain stays under 1/(2 l_R) for every r.
 */

#include "RmpccAssociatedPhase.h"
#include "RmpccPhaseResidual.h"

#include <Math/MathFunctions.h>

#include <Eigen/Core>
#include <Eigen/SVD>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>

namespace {

using RobotLibrary::Control::RmpccContourResidualGeometry;
using RobotLibrary::Control::RmpccPhaseAssociation;
using RobotLibrary::Control::RmpccPoseArcTable;
using RobotLibrary::Control::RmpccStateVector;

constexpr double kRegularization = 1e-10;
constexpr double kTolerance = 1e-12;
constexpr double kLinearTangent = 0.176;
constexpr double kFiniteDifferenceStep = 1e-7;

double gRatio = 1.0;
double gRotationLength = 0.20;
Eigen::Vector3d gRotationAxis = Eigen::Vector3d::UnitZ();

Eigen::Vector<double,6> screw()
{
    Eigen::Vector<double,6> value;
    value.head<3>() = kLinearTangent * Eigen::Vector3d::UnitX();
    value.tail<3>() = gRatio * kLinearTangent * gRotationAxis;
    return value;
}

Eigen::Matrix4d reference(const double progress)
{
    return RobotLibrary::Math::se3_exponential(screw() * progress);
}

Eigen::Vector<double,6> tangent(const double) { return screw(); }

Eigen::Matrix<double,6,6> metric()
{
    Eigen::Matrix<double,6,6> value = Eigen::Matrix<double,6,6>::Identity();
    value.diagonal().tail<3>().setConstant(gRotationLength * gRotationLength);
    return value;
}

const char *name(const RmpccPhaseAssociation association)
{
    switch(association)
    {
        case RmpccPhaseAssociation::MetricScrew: return "MetricScrew";
        case RmpccPhaseAssociation::TaskPointXYZ: return "TaskPointXYZ";
        case RmpccPhaseAssociation::TaskPoseFeature: return "TaskPoseFeature";
    }
    return "Unknown";
}

/// Spectral norm of d(contour rotation)/d(position error), in rad/m.
double associated_coupling(const RmpccStateVector &state,
                           const RmpccPhaseAssociation association,
                           const RmpccContourResidualGeometry geometry)
{
    RmpccPoseArcTable poseArc;
    poseArc.build(metric(), tangent, 4097);
    const auto linearized =
        RobotLibrary::Control::rmpcc_linearize_associated_residuals(
            state, geometry, metric(), association, kRegularization,
            kTolerance, gRotationLength, kFiniteDifferenceStep, poseArc,
            reference, tangent);
    return Eigen::JacobiSVD<Eigen::MatrixXd>(
        Eigen::MatrixXd(linearized.contourJacobian.block<3,3>(3,0)))
        .singularValues()(0);
}

double associated_phase_correction(const RmpccStateVector &state,
                                   const RmpccPhaseAssociation association)
{
    RmpccPoseArcTable poseArc;
    poseArc.build(metric(), tangent, 4097);
    return RobotLibrary::Control::rmpcc_associated_residuals(
               state,
               RmpccContourResidualGeometry::AssociatedDecoupledCartesianSO3,
               metric(), association, kRegularization, kTolerance,
               gRotationLength, poseArc, reference, tangent)
        .context.phaseCorrection;
}

double local_phase_correction(const RmpccStateVector &state,
                              const RmpccPhaseAssociation association)
{
    return RobotLibrary::Control::rmpcc_phase_residuals(
               state, metric(), association, kRegularization, kTolerance,
               gRotationLength, reference, tangent)
        .phaseCorrection;
}

} // namespace

int main()
{
    constexpr std::array<double,5> ratios{{1.0, 5.0, 10.0, 20.0, 50.0}};
    constexpr std::array<double,3> rotationLengths{{0.05, 0.20, 0.40}};
    const std::array<Eigen::Vector3d,3> axes{{
        Eigen::Vector3d::UnitX(),                                                                    // q(a) = 1, the least observable case
        Eigen::Vector3d::UnitZ(),                                                                    // q(a) = 2
        Eigen::Vector3d(1.0, 2.0, 3.0).normalized()}};
    constexpr std::array<RmpccPhaseAssociation,3> associations{{
        RmpccPhaseAssociation::MetricScrew,
        RmpccPhaseAssociation::TaskPointXYZ,
        RmpccPhaseAssociation::TaskPoseFeature}};

    // A state off the path, at the closed-loop working point: the conditioning
    // claim has to hold where the controller lives, not only at eta = 0.
    RmpccStateVector working = RmpccStateVector::Zero();
    working(0) = 6e-4;                                                                               // position error [m]
    working(5) = 5e-3;                                                                               // rotation error [rad]
    working(6) = 0.5;
    RmpccStateVector onPath = RmpccStateVector::Zero();
    onPath(6) = 0.5;

    int failures = 0;
    double worstFeatureCoupling = 0.0;
    double worstDisagreement = 0.0;

    for(const double rotationLength : rotationLengths)
    {
        for(const Eigen::Vector3d &axis : axes)
        {
            for(const double ratio : ratios)
            {
                gRatio = ratio;
                gRotationAxis = axis;
                gRotationLength = rotationLength;

                // 1. One implementation: both families produce the same delta_s.
                for(const auto association : associations)
                {
                    for(const RmpccStateVector &state : {onPath, working})
                    {
                        const double localValue =
                            local_phase_correction(state, association);
                        const double associatedValue =
                            associated_phase_correction(state, association);
                        const double disagreement =
                            std::abs(localValue - associatedValue);
                        worstDisagreement =
                            std::max(worstDisagreement, disagreement);
                        if(disagreement > 1e-15)
                        {
                            std::printf(
                                "local/associated delta_s disagree: %s r=%g l_R=%g "
                                "local=%.17g associated=%.17g\n",
                                name(association), ratio, rotationLength,
                                localValue, associatedValue);
                            ++failures;
                        }
                    }
                }

                for(const auto geometry :
                    {RmpccContourResidualGeometry::ExactAssociatedUnifiedSE3,
                     RmpccContourResidualGeometry::AssociatedDecoupledCartesianSO3})
                {
                    // 2. Known-bad baseline, characterised where it is exact:
                    // on the path the position-to-rotation gain is precisely
                    // the tangent ratio. Off the path it drifts by O(|eta|),
                    // so the equality is asserted at eta = 0.
                    const double taskPoint = associated_coupling(
                        onPath, RmpccPhaseAssociation::TaskPointXYZ, geometry);
                    if(std::abs(taskPoint - ratio) / ratio > 1e-6)
                    {
                        std::printf(
                            "associated TaskPointXYZ gain mismatch: r=%g l_R=%g "
                            "measured=%.10g expected=%g\n",
                            ratio, rotationLength, taskPoint, ratio);
                        ++failures;
                    }

                    // 3. Repaired mode: the bound must hold everywhere the
                    // controller operates, not only at the linearisation point.
                    const double bound = 1.0 / (2.0 * rotationLength);
                    for(const RmpccStateVector &state : {onPath, working})
                    {
                        const double feature = associated_coupling(
                            state, RmpccPhaseAssociation::TaskPoseFeature,
                            geometry);
                        worstFeatureCoupling =
                            std::max(worstFeatureCoupling, feature / bound);
                        if(not std::isfinite(feature) or feature > bound + 1e-6)
                        {
                            std::printf(
                                "associated TaskPoseFeature exceeded bound: r=%g "
                                "l_R=%g measured=%.10g bound=%.10g\n",
                                ratio, rotationLength, feature, bound);
                            ++failures;
                        }
                    }
                }
            }
        }
    }

    std::printf("worst local/associated delta_s disagreement = %.3g\n",
                worstDisagreement);
    std::printf("worst TaskPoseFeature coupling / (1/(2 l_R)) = %.6g\n",
                worstFeatureCoupling);
    return failures == 0 ? 0 : 1;
}
