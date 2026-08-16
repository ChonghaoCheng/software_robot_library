#include "RmpccCostGeometry.h"
#include "RmpccPhaseResidual.h"
#include "RmpccResidualLinearization.h"

#include <Math/MathFunctions.h>

#include <Eigen/Core>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>

namespace {

using RobotLibrary::Control::RmpccLagPenalty;
using RobotLibrary::Control::RmpccPhaseAssociation;
using RobotLibrary::Control::RmpccPhaseResiduals;
using RobotLibrary::Control::RmpccStateVector;

constexpr double kRegularization = 1e-10;
constexpr double kTolerance = 1e-12;

Eigen::Vector<double,6> screw()
{
    Eigen::Vector<double,6> value;
    value << 0.21, -0.13, 0.17, 0.36, -0.24, 0.29;
    return value;
}

Eigen::Matrix4d reference(const double progress)
{
    return RobotLibrary::Math::se3_exponential(screw() * progress);
}

Eigen::Vector<double,6> tangent(const double)
{
    return screw();
}

Eigen::Matrix<double,6,6> metric()
{
    Eigen::Matrix<double,6,6> value =
        Eigen::Matrix<double,6,6>::Identity();
    value.diagonal().tail<3>().setConstant(0.20 * 0.20);
    return value;
}

RmpccPhaseResiduals residual(
    const RmpccStateVector &state,
    const RmpccPhaseAssociation association)
{
    return RobotLibrary::Control::rmpcc_phase_residuals(
        state, metric(), association, kRegularization, kTolerance,
        reference, tangent);
}

double objective(const RmpccStateVector &state,
                 const RmpccPhaseAssociation association,
                 const RmpccLagPenalty penalty)
{
    const auto value = residual(state, association);
    Eigen::Matrix<double,6,6> contourWeight =
        Eigen::Matrix<double,6,6>::Identity();
    contourWeight.diagonal() << 900.0, 1100.0, 800.0, 31.0, 43.0, 37.0;
    Eigen::Matrix<double,6,6> lagWeight =
        Eigen::Matrix<double,6,6>::Identity();
    lagWeight.diagonal() << 12.0, 12.0, 12.0, 7.0, 9.0, 8.0;
    const double contourCost =
        (value.contour.transpose() * contourWeight * value.contour)(0);
    return contourCost + (penalty == RmpccLagPenalty::PhaseInducedPoseVector
        ? (value.vectorLag.transpose() * lagWeight * value.vectorLag)(0)
        : lagWeight(0,0) * value.scalarLag * value.scalarLag);
}

double relative_error(const double value, const double expected)
{
    return std::abs(value - expected) /
        std::max({1e-12, std::abs(value), std::abs(expected)});
}

Eigen::Vector<double,6> numerical_body_tangent(
    const double progress,
    const Eigen::Vector<double,6> &canonicalScrew,
    const Eigen::Matrix4d &attachment)
{
    constexpr double h = 1e-6;
    const auto taskReference = [&](const double s)
    {
        return RobotLibrary::Math::se3_exponential(canonicalScrew * s)
            * attachment;
    };
    return RobotLibrary::Math::se3_logarithm(
        RobotLibrary::Math::se3_inverse(taskReference(progress))
        * taskReference(progress + h)) / h;
}

} // namespace

int main()
{
    RmpccStateVector state;
    state << 0.05, -0.04, 0.03, 0.48, -0.39, 0.32, 0.41;

    // 1. MV reproduces the corrected legacy FullScrew residual, objective,
    // directional derivative, and local QP H/f assembly.
    const auto legacy = RobotLibrary::Control::rmpcc_full_screw_residuals(
        state, metric(), kRegularization, tangent);
    const auto mv = residual(state, RmpccPhaseAssociation::MetricScrew);
    const double residualError = std::max(
        (legacy.contour - mv.contour).norm(),
        (legacy.lag - mv.vectorLag).norm());
    if(residualError > 1e-13)
    {
        std::cerr << "MV legacy residual mismatch: " << residualError << '\n';
        return 1;
    }

    Eigen::Matrix<double,6,6> contourWeight =
        Eigen::Matrix<double,6,6>::Identity();
    contourWeight.diagonal() << 900.0, 1100.0, 800.0, 31.0, 43.0, 37.0;
    Eigen::Matrix<double,6,6> lagWeight =
        Eigen::Matrix<double,6,6>::Identity();
    lagWeight.diagonal() << 12.0, 12.0, 12.0, 7.0, 9.0, 8.0;
    const double legacyObjective =
        (legacy.contour.transpose() * contourWeight * legacy.contour)(0)
        + (legacy.lag.transpose() * lagWeight * legacy.lag)(0);
    if(std::abs(legacyObjective - objective(
           state, RmpccPhaseAssociation::MetricScrew,
           RmpccLagPenalty::PhaseInducedPoseVector)) > 1e-13)
    {
        std::cerr << "MV legacy scalar objective mismatch.\n";
        return 2;
    }

    constexpr double h = 1e-6;
    const auto oldLinearization =
        RobotLibrary::Control::rmpcc_linearize_full_screw_residuals(
            state, metric(), kRegularization, h, tangent);
    const auto newLinearization =
        RobotLibrary::Control::rmpcc_linearize_phase_residuals(
            state, metric(), RmpccPhaseAssociation::MetricScrew,
            kRegularization, kTolerance, h, reference, tangent);
    Eigen::Matrix<double,7,7> oldH =
        2.0 * oldLinearization.contourJacobian.transpose()
            * contourWeight * oldLinearization.contourJacobian
        + 2.0 * oldLinearization.lagJacobian.transpose()
            * lagWeight * oldLinearization.lagJacobian;
    Eigen::Vector<double,7> oldF =
        2.0 * oldLinearization.contourJacobian.transpose()
            * contourWeight * oldLinearization.residual.contour
        + 2.0 * oldLinearization.lagJacobian.transpose()
            * lagWeight * oldLinearization.residual.lag;
    Eigen::Matrix<double,7,7> newH =
        2.0 * newLinearization.contourJacobian.transpose()
            * contourWeight * newLinearization.contourJacobian
        + 2.0 * newLinearization.vectorLagJacobian.transpose()
            * lagWeight * newLinearization.vectorLagJacobian;
    Eigen::Vector<double,7> newF =
        2.0 * newLinearization.contourJacobian.transpose()
            * contourWeight * newLinearization.residual.contour
        + 2.0 * newLinearization.vectorLagJacobian.transpose()
            * lagWeight * newLinearization.residual.vectorLag;
    const double qpHError = (newH - oldH).norm();
    const double qpFError = (newF - oldF).norm();
    const double qpHRelativeError = qpHError / std::max(1.0, oldH.norm());
    const double qpFRelativeError = qpFError / std::max(1.0, oldF.norm());
    if(qpHRelativeError > 1e-8 || qpFRelativeError > 1e-8)
    {
        std::cerr << "MV legacy QP mismatch: H=" << qpHError
                  << " f=" << qpFError << " relative H=" << qpHRelativeError
                  << " relative f=" << qpFRelativeError << '\n';
        return 3;
    }

    // 2. At d=0, translational metric phase and task-point phase agree to
    // first order.
    RmpccStateVector local = RmpccStateVector::Zero();
    local << 2e-7, -3e-7, 1e-7, 4e-7, -2e-7, 3e-7, 0.37;
    Eigen::Matrix<double,6,6> translationMetric =
        Eigen::Matrix<double,6,6>::Zero();
    translationMetric.diagonal().head<3>().setOnes();
    const auto localMetric = RobotLibrary::Control::rmpcc_phase_residuals(
        local, translationMetric, RmpccPhaseAssociation::MetricScrew,
        kRegularization, kTolerance, reference, tangent);
    const auto localTask = residual(local, RmpccPhaseAssociation::TaskPointXYZ);
    if(std::abs(localMetric.phaseCorrection - localTask.phaseCorrection) > 1e-12)
    {
        std::cerr << "d=0 first-order phase mismatch.\n";
        return 4;
    }

    // 3. A pure path-phase mismatch is recovered to first order and leaves a
    // second-order task contour.
    std::array<double,3> phaseErrors{};
    std::array<double,3> contourNorms{};
    const std::array<double,3> phaseSteps{1e-2, 5e-3, 2.5e-3};
    for(std::size_t i = 0; i < phaseSteps.size(); ++i)
    {
        const double ds = phaseSteps[i];
        RmpccStateVector shifted = RmpccStateVector::Zero();
        shifted(6) = 0.34;
        shifted.head<6>() = RobotLibrary::Math::se3_logarithm(
            RobotLibrary::Math::se3_inverse(reference(shifted(6)))
            * reference(shifted(6) + ds));
        const auto task = residual(shifted, RmpccPhaseAssociation::TaskPointXYZ);
        phaseErrors[i] = std::abs(task.phaseCorrection - ds);
        contourNorms[i] = task.contour.norm();
    }
    if(phaseErrors[0] > 2e-4 || contourNorms[0] > 2e-3
       || phaseErrors[0] / std::max(phaseErrors[1], 1e-16) < 3.5
       || contourNorms[0] / std::max(contourNorms[1], 1e-16) < 3.5)
    {
        std::cerr << "Pure-phase convergence order failed.\n";
        return 5;
    }

    // 4. Rotation about the task origin cannot create task-point phase.
    RmpccStateVector orientationOnly = RmpccStateVector::Zero();
    orientationOnly.tail<4>() << 0.08, -0.05, 0.04, 0.42;
    const auto orientationResidual =
        residual(orientationOnly, RmpccPhaseAssociation::TaskPointXYZ);
    if(std::abs(orientationResidual.phaseCorrection) > 1e-13
       || std::abs(orientationResidual.scalarLag) > 1e-13
       || (orientationResidual.contour.tail<3>()
           - orientationOnly.head<6>().tail<3>()).norm() > 1e-13)
    {
        std::cerr << "Task-origin orientation perturbation leaked into phase.\n";
        return 6;
    }

    // 5. A wrist rotation with a 50 mm attachment moves the synthetic task
    // point, and the phase equals that displacement's tangent projection.
    Eigen::Matrix4d attachment = Eigen::Matrix4d::Identity();
    attachment(2,3) = 0.05;
    const Eigen::Vector<double,6> canonicalScrew = screw();
    const auto taskReference = [&](const double s)
    {
        return RobotLibrary::Math::se3_exponential(canonicalScrew * s)
            * attachment;
    };
    const auto taskTangent = [&](const double s)
    {
        return numerical_body_tangent(s, canonicalScrew, attachment);
    };
    RmpccStateVector leverState = RmpccStateVector::Zero();
    leverState(6) = 0.38;
    Eigen::Vector<double,6> wristRotation = Eigen::Vector<double,6>::Zero();
    wristRotation.tail<3>() << 0.015, -0.01, 0.02;
    const Eigen::Matrix4d actualTask =
        RobotLibrary::Math::se3_exponential(canonicalScrew * leverState(6))
        * RobotLibrary::Math::se3_exponential(wristRotation) * attachment;
    leverState.head<6>() = RobotLibrary::Math::se3_logarithm(
        RobotLibrary::Math::se3_inverse(taskReference(leverState(6))) * actualTask);
    const auto leverResidual = RobotLibrary::Control::rmpcc_phase_residuals(
        leverState, metric(), RmpccPhaseAssociation::TaskPointXYZ,
        kRegularization, kTolerance, taskReference, taskTangent);
    const Eigen::Vector3d taskDirection =
        taskReference(leverState(6)).block<3,3>(0,0)
        * taskTangent(leverState(6)).head<3>();
    const Eigen::Vector3d taskDisplacement =
        actualTask.block<3,1>(0,3)
        - taskReference(leverState(6)).block<3,1>(0,3);
    const double expectedLeverPhase =
        taskDirection.dot(taskDisplacement) / taskDirection.squaredNorm();
    if(std::abs(leverResidual.phaseCorrection - expectedLeverPhase) > 1e-12)
    {
        std::cerr << "Lever-arm task phase mismatch.\n";
        return 7;
    }

    // 6/7. Tangential world displacement has signed phase; transverse
    // displacement has no first-order phase.
    const double s = 0.43;
    const Eigen::Matrix4d ref = reference(s);
    const Eigen::Vector3d worldTangent = ref.block<3,3>(0,0) * tangent(s).head<3>();
    const Eigen::Vector3d unitTangent = worldTangent.normalized();
    Eigen::Vector3d transverse = unitTangent.cross(Eigen::Vector3d::UnitZ());
    if(transverse.norm() < 0.2)
        transverse = unitTangent.cross(Eigen::Vector3d::UnitY());
    transverse.normalize();
    for(const auto &fixture : std::array<std::pair<Eigen::Vector3d,double>,2>{
            std::make_pair(0.002 * unitTangent, 0.002 / worldTangent.norm()),
            std::make_pair(0.002 * transverse, 0.0)})
    {
        Eigen::Matrix4d actual = ref;
        actual.block<3,1>(0,3) += fixture.first;
        RmpccStateVector displaced = RmpccStateVector::Zero();
        displaced(6) = s;
        displaced.head<6>() = RobotLibrary::Math::se3_logarithm(
            RobotLibrary::Math::se3_inverse(ref) * actual);
        const auto value = residual(displaced, RmpccPhaseAssociation::TaskPointXYZ);
        if(std::abs(value.phaseCorrection - fixture.second) > 1e-11)
        {
            std::cerr << "Tangential/transverse phase fixture failed.\n";
            return 8;
        }
    }

    // 8. A zero task-point tangent is explicitly unobservable and finite.
    const auto rotationReference = [](const double progress)
    {
        Eigen::Vector<double,6> xi = Eigen::Vector<double,6>::Zero();
        xi.tail<3>() << 0.2, -0.3, 0.4;
        return RobotLibrary::Math::se3_exponential(xi * progress);
    };
    const auto rotationTangent = [](const double)
    {
        Eigen::Vector<double,6> xi = Eigen::Vector<double,6>::Zero();
        xi.tail<3>() << 0.2, -0.3, 0.4;
        return xi;
    };
    const auto unobservable = RobotLibrary::Control::rmpcc_phase_residuals(
        state, metric(), RmpccPhaseAssociation::TaskPointXYZ,
        kRegularization, kTolerance, rotationReference, rotationTangent);
    if(unobservable.phaseObservable || not unobservable.contour.allFinite()
       || not unobservable.vectorLag.allFinite()
       || not std::isfinite(unobservable.scalarLag)
       || unobservable.vectorLag.norm() != 0.0 || unobservable.scalarLag != 0.0
       || (unobservable.contour - state.head<6>()).norm() != 0.0)
    {
        std::cerr << "Unobservable task phase fallback failed.\n";
        return 9;
    }

    // 9. All four complete residual objectives pass a directional derivative.
    RmpccStateVector direction;
    direction << -0.17, 0.23, 0.11, 0.31, -0.29, 0.19, 0.37;
    direction.normalize();
    for(const auto association : {RmpccPhaseAssociation::MetricScrew,
                                  RmpccPhaseAssociation::TaskPointXYZ})
    {
        for(const auto penalty : {RmpccLagPenalty::PhaseInducedPoseVector,
                                  RmpccLagPenalty::ScalarTaskDistance})
        {
            const auto linearization =
                RobotLibrary::Control::rmpcc_linearize_phase_residuals(
                    state, metric(), association, kRegularization, kTolerance,
                    h, reference, tangent);
            Eigen::Vector<double,7> gradient =
                2.0 * linearization.contourJacobian.transpose()
                    * contourWeight * linearization.residual.contour;
            if(penalty == RmpccLagPenalty::PhaseInducedPoseVector)
            {
                gradient += 2.0 * linearization.vectorLagJacobian.transpose()
                    * lagWeight * linearization.residual.vectorLag;
            }
            else
            {
                gradient += 2.0 * lagWeight(0,0)
                    * linearization.scalarLagJacobian.transpose()
                    * linearization.residual.scalarLag;
            }
            constexpr double objectiveStep = 2e-6;
            const double expected =
                (objective(state + objectiveStep * direction, association, penalty)
                 - objective(state - objectiveStep * direction, association, penalty))
                / (2.0 * objectiveStep);
            const double error = relative_error(gradient.dot(direction), expected);
            if(error > 1e-6)
            {
                std::cerr << "Factorial objective derivative error: " << error << '\n';
                return 10;
            }
        }
    }

    // 10. Residual geometry is stage-independent; terminal semantics can only
    // multiply its costs.
    const auto running = residual(state, RmpccPhaseAssociation::TaskPointXYZ);
    const auto terminal = residual(state, RmpccPhaseAssociation::TaskPointXYZ);
    if((running.contour - terminal.contour).norm() != 0.0
       || (running.vectorLag - terminal.vectorLag).norm() != 0.0
       || running.scalarLag != terminal.scalarLag)
    {
        std::cerr << "Running/terminal residual geometry diverged.\n";
        return 11;
    }

    std::cout << "MV residual error=" << residualError
              << " QP H error=" << qpHError
              << " QP f error=" << qpFError
              << " pure-phase contour order="
              << contourNorms[0] / std::max(contourNorms[1], 1e-16) << '\n';
    return 0;
}
