/**
 * @file    SerialLinkTimeIndexedMPC.cpp
 * @brief   Behaviour-neutral diagnostics recording shared by time-indexed MPCs.
 */

#include <Control/TrajectoryTracking/SerialLinkTimeIndexedMPC.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace RobotLibrary { namespace Control {

void
SerialLinkTimeIndexedMPC::record_time_indexed_diagnostics(
    const Eigen::Vector<double,6> &commandedTwistInBase,
    const Eigen::Vector<double,6> &clampFrameTwist,
    const Eigen::VectorXd &jointCommand,
    const double referenceTime,
    const double controllerDt,
    const double maxLinearSpeed,
    const double maxAngularSpeed,
    const double solveTimeSeconds)
{
    TimeIndexedMpcDiagnostics &d = _timeIndexedDiagnostics;

    d.invocationIndex += 1.0;
    d.controllerDt = controllerDt;
    d.referenceTime = referenceTime;
    d.qpStatus = 1.0;                 // A returned invocation solved; failures throw.
    d.commandedTwist = commandedTwistInBase;
    d.clampFrameTwist = clampFrameTwist;
    d.effectiveLinearLimit = maxLinearSpeed;
    d.effectiveAngularLimit = maxAngularSpeed;
    d.solveTimeSeconds = solveTimeSeconds;

    // Realized twist from the same Jacobian and the joint command actually
    // returned, so the error is the incurred one rather than a reconstruction.
    if(_jacobianMatrix.cols() == jointCommand.size())
    {
        d.realizedTwist = _jacobianMatrix * jointCommand;
    }
    else
    {
        d.realizedTwist.setZero();
    }
    d.twistRealizationError = (d.realizedTwist - d.commandedTwist).norm();

    // A Cartesian axis is saturated when the post-clamp command sits on the
    // clamp, compared in the axes the clamp was actually applied in. The clamp
    // is applied by the caller, so equality is exact.
    d.linearLimitActive = false;
    d.angularLimitActive = false;
    const double clampTolerance = 1e-12;
    for(int i = 0; i < 3; ++i)
    {
        if(std::abs(std::abs(clampFrameTwist(i)) - maxLinearSpeed) <= clampTolerance)
            d.linearLimitActive = true;
        if(std::abs(std::abs(clampFrameTwist(i + 3)) - maxAngularSpeed) <= clampTolerance)
            d.angularLimitActive = true;
    }

    // A joint is saturated when its command sits on the same control limit the
    // resolved-rate QP was constrained by.
    // The margin is recorded alongside the flag so the flag can be re-derived
    // from published data, without exposing the protected limits themselves.
    const double jointTolerance = 1e-9;
    double margin = std::numeric_limits<double>::infinity();
    for(Eigen::Index i = 0; i < jointCommand.size(); ++i)
    {
        const RobotLibrary::Model::Limits limits =
            compute_control_limits(static_cast<unsigned int>(i));
        margin = std::min(margin, jointCommand(i) - limits.lower);
        margin = std::min(margin, limits.upper - jointCommand(i));
    }
    if(not std::isfinite(margin)) margin = 0.0;
    d.jointLimitMargin = margin;
    d.jointLimitActive = margin <= jointTolerance;

    d.finite = d.commandedTwist.allFinite() and d.clampFrameTwist.allFinite()
               and d.realizedTwist.allFinite()
               and std::isfinite(d.twistRealizationError)
               and std::isfinite(d.controllerDt) and std::isfinite(d.referenceTime)
               and std::isfinite(d.jointLimitMargin)
               and jointCommand.allFinite();
}

} } // namespace RobotLibrary::Control
