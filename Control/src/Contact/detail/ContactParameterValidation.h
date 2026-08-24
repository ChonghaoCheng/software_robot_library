#ifndef CONTROL_CONTACT_DETAIL_CONTACT_PARAMETER_VALIDATION_H
#define CONTROL_CONTACT_DETAIL_CONTACT_PARAMETER_VALIDATION_H

#include <Control/Contact/ContactDataStructures.h>

#include <cmath>
#include <stdexcept>

namespace RobotLibrary { namespace Control { namespace detail {

inline bool finite_nonnegative(const double value)
{
    return std::isfinite(value) && value >= 0.0;
}

inline void validate_contact_parameters(const ContactParameters &parameters)
{
    if(!std::isfinite(parameters.forceResponseGain) || parameters.forceResponseGain <= 0.0)
        throw std::invalid_argument("ContactParameters.forceResponseGain must be finite and positive.");
    if(!finite_nonnegative(parameters.targetForce))
        throw std::invalid_argument("ContactParameters.targetForce must be finite and nonnegative.");
    if(!finite_nonnegative(parameters.forceTolerance))
        throw std::invalid_argument("ContactParameters.forceTolerance must be finite and nonnegative.");
    if(!finite_nonnegative(parameters.forceWeight))
        throw std::invalid_argument("ContactParameters.forceWeight must be finite and nonnegative.");
    if(!finite_nonnegative(parameters.slackWeight))
        throw std::invalid_argument("ContactParameters.slackWeight must be finite and nonnegative.");
    if(!finite_nonnegative(parameters.maxForceSlack))
        throw std::invalid_argument("ContactParameters.maxForceSlack must be finite and nonnegative.");
    if(!finite_nonnegative(parameters.tangentPositionWeight))
        throw std::invalid_argument("ContactParameters.tangentPositionWeight must be finite and nonnegative.");
    if(!finite_nonnegative(parameters.pathLagPositionWeight))
        throw std::invalid_argument("ContactParameters.pathLagPositionWeight must be finite and nonnegative.");
    if(!finite_nonnegative(parameters.normalPositionWeight))
        throw std::invalid_argument("ContactParameters.normalPositionWeight must be finite and nonnegative.");
    if(!finite_nonnegative(parameters.orientationWeight))
        throw std::invalid_argument("ContactParameters.orientationWeight must be finite and nonnegative.");
    if(!finite_nonnegative(parameters.velocityWeight))
        throw std::invalid_argument("ContactParameters.velocityWeight must be finite and nonnegative.");
    if(!finite_nonnegative(parameters.deltaVelocityWeight))
        throw std::invalid_argument("ContactParameters.deltaVelocityWeight must be finite and nonnegative.");
    if(!finite_nonnegative(parameters.maxBoardAcceleration))
        throw std::invalid_argument("ContactParameters.maxBoardAcceleration must be finite and nonnegative.");
    if(!parameters.normalAxisInBoard.allFinite() || parameters.normalAxisInBoard.norm() <= 1e-9)
        throw std::invalid_argument("ContactParameters.normalAxisInBoard must be finite and nonzero.");
    if(!parameters.tangentAxisInBoard.allFinite() || parameters.tangentAxisInBoard.norm() <= 1e-9)
        throw std::invalid_argument("ContactParameters.tangentAxisInBoard must be finite and nonzero.");
    const Eigen::Vector3d normal = parameters.normalAxisInBoard.normalized();
    const Eigen::Vector3d tangent = parameters.tangentAxisInBoard.normalized();
    if(std::abs(normal.dot(tangent)) > 1.0 - 1e-6)
        throw std::invalid_argument("ContactParameters normal and tangent axes must not be parallel.");

    switch(parameters.mode)
    {
        case ContactMode::Disabled:
        case ContactMode::Loss:
        case ContactMode::Constraint:
            break;
        default:
            throw std::invalid_argument("ContactParameters.mode is invalid.");
    }
}

inline double validate_measured_normal_force(const double force)
{
    if(!finite_nonnegative(force))
        throw std::invalid_argument("Measured normal force must be finite and nonnegative.");
    return force;
}

} } } // namespace RobotLibrary::Control::detail

#endif
