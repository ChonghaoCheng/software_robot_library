/**
 * @file    AdmittanceContactController.cpp
 * @brief   Baseline admittance-style normal-force contact controller.
 */

#include <Control/Contact/AdmittanceContactController.h>

#include <algorithm>
#include <cmath>

namespace RobotLibrary { namespace Control {

double
AdmittanceContactController::apply_deadband(double value, double deadband)
{
    const double width = std::max(0.0, deadband);
    if(std::abs(value) <= width)
    {
        return 0.0;
    }

    return std::copysign(std::abs(value) - width, value);
}

ContactCommand
AdmittanceContactController::compute(const ContactState &state,
                                     const ContactReference &reference,
                                     double dt)
{
    (void)dt;

    ContactCommand command;
    _diagnostics = ContactDiagnostics{};

    const bool normalValid = normal_is_valid(state.normalBase);
    _diagnostics.normalValid = normalValid;
    _diagnostics.measuredNormalForce = std::max(0.0, state.measuredNormalForce);
    _diagnostics.targetNormalForce = std::max(0.0, reference.targetNormalForce);

    if(!reference.enabled || !normalValid)
    {
        return command;
    }

    const Eigen::Vector3d normal = unit_normal(state.normalBase);
    const double maxNormalVelocity = std::max(0.0, reference.maxNormalVelocity);
    const double forceResponseGain = std::max(1e-9, reference.forceResponseGain);
    const double surfaceNormalVelocity = normal.dot(state.surfaceVelocityBase);
    const double forceError = _diagnostics.targetNormalForce - _diagnostics.measuredNormalForce;

    command.forceError = forceError;
    command.active = state.inContact || reference.approachWhenNotInContact;
    command.approachActive =
        !state.inContact && reference.approachWhenNotInContact && _diagnostics.targetNormalForce > 0.0;

    double relativeNormalVelocity = 0.0;
    if(command.approachActive)
    {
        relativeNormalVelocity = std::abs(reference.approachVelocity);
    }
    else if(command.active)
    {
        const double controlledForceError = apply_deadband(forceError, reference.forceDeadband);
        relativeNormalVelocity =
            reference.forceVelocityGain * controlledForceError / forceResponseGain;
    }

    relativeNormalVelocity =
        std::clamp(relativeNormalVelocity, -maxNormalVelocity, maxNormalVelocity);

    command.desiredRelativeNormalVelocity = relativeNormalVelocity;
    command.desiredRobotNormalVelocity = surfaceNormalVelocity + relativeNormalVelocity;
    command.normalVelocityBase = normal * command.desiredRobotNormalVelocity;
    command.releaseActive = command.active && relativeNormalVelocity < 0.0;

    _diagnostics.forceError = command.forceError;
    _diagnostics.surfaceNormalVelocity = surfaceNormalVelocity;
    _diagnostics.desiredRelativeNormalVelocity = command.desiredRelativeNormalVelocity;
    _diagnostics.desiredRobotNormalVelocity = command.desiredRobotNormalVelocity;
    _diagnostics.active = command.active;
    _diagnostics.approachActive = command.approachActive;
    _diagnostics.releaseActive = command.releaseActive;

    return command;
}

} } // namespace
