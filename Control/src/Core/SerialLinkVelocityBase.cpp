/**
 * @file    SerialLinkVelocityBase.cpp
 * @brief   Shared resolved-rate control for serial link velocity controllers.
 */

#include <Control/Core/SerialLinkVelocityBase.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace RobotLibrary { namespace Control {

SerialLinkVelocityBase::SerialLinkVelocityBase(
    std::shared_ptr<RobotLibrary::Model::KinematicTree> model,
    const std::string &endpointName,
    const RobotLibrary::Control::SerialLinkParameters &parameters)
: SerialLinkBase(model, endpointName, parameters)
{
}

Eigen::VectorXd
SerialLinkVelocityBase::resolve_endpoint_twist(const Eigen::Vector<double,6> &twist)
{
    return resolve_endpoint_motion(twist);
}

Eigen::VectorXd
SerialLinkVelocityBase::resolve_endpoint_motion(const Eigen::Vector<double,6> &endpointMotion)
{
    using namespace Eigen;

    const unsigned int numJoints = _model->number_of_joints();
    VectorXd startPoint = _model->joint_velocities();
    VectorXd lowerBound(numJoints), upperBound(numJoints);

    for (unsigned int i = 0; i < numJoints; ++i)
    {
        const auto &[lower, upper] = compute_control_limits(i);
        lowerBound[i] = lower;
        upperBound[i] = upper;

        startPoint[i] = std::clamp(startPoint[i], lower + 1e-03, upper - 1e-03);
    }

    const VectorXd manipulabilityGradient = manipulability_gradient();
    _constraintMatrix.row(2 * numJoints) = -manipulabilityGradient.transpose();
    _constraintVector.head(numJoints) = upperBound;
    _constraintVector.segment(numJoints, numJoints) = -lowerBound;
    _constraintVector(2 * numJoints) =
        (_manipulability - _minManipulability) * 100 * std::sqrt(_controlFrequency);

    VectorXd controlVelocity = VectorXd::Zero(numJoints);

    if (not is_singular())
    {
        if (numJoints <= 6)
        {
            controlVelocity = QPSolver<double>::solve(
                _jacobianMatrix.transpose() * _jacobianMatrix,
               -_jacobianMatrix.transpose() * endpointMotion,
                _constraintMatrix,
                _constraintVector,
                startPoint);
        }
        else
        {
            if (not _redundantTaskSet)
            {
                _redundantTask =
                    manipulabilityGradient * std::sqrt(_controlFrequency) / 10.0;
                _redundantTaskSet = false;
            }

            controlVelocity = QPSolver<double>::constrained_least_squares(
                _redundantTask,
                _model->joint_inertia_matrix(),
                _jacobianMatrix,
                endpointMotion,
                _constraintMatrix,
                _constraintVector,
                startPoint);
        }
    }
    else
    {
        const double dampingFactor =
            std::pow(1.0 - _manipulability / _minManipulability, 2.0) * 0.01;

        MatrixXd H = _jacobianMatrix.transpose() * _jacobianMatrix;
        H.diagonal().array() += dampingFactor;

        controlVelocity = QPSolver<double>::solve(
            H,
           -_jacobianMatrix.transpose() * endpointMotion,
            _constraintMatrix.block(0, 0, 2 * numJoints, numJoints),
            _constraintVector.head(2 * numJoints),
            startPoint);
    }

    return controlVelocity;
}

Eigen::VectorXd
SerialLinkVelocityBase::track_joint_trajectory(
    const Eigen::VectorXd &desiredPosition,
    const Eigen::VectorXd &desiredVelocity,
    const Eigen::VectorXd &desiredAcceleration)
{
    (void)desiredAcceleration;

    const unsigned int numJoints = _model->number_of_joints();

    if (desiredPosition.size() != numJoints or desiredVelocity.size() != numJoints)
    {
        throw std::invalid_argument(
            "[ERROR] [SERIAL LINK VELOCITY] track_joint_trajectory(): "
            "Incorrect size for input arguments. This robot has "
            + std::to_string(numJoints) + " joints, but "
            "the position argument had " + std::to_string(desiredPosition.size())
            + " elements, and the velocity argument had "
            + std::to_string(desiredVelocity.size()) + " elements.");
    }

    Eigen::VectorXd velocityControl(numJoints);

    for (unsigned int i = 0; i < numJoints; ++i)
    {
        velocityControl(i) =
            desiredVelocity(i)
            + _jointPositionGains[i]
                * (desiredPosition(i) - _model->joint_positions()[i]);

        const RobotLibrary::Model::Limits controlLimits = compute_control_limits(i);

        if (velocityControl(i) <= controlLimits.lower)
        {
            velocityControl(i) = controlLimits.lower + 1e-03;
        }
        else if (velocityControl(i) >= controlLimits.upper)
        {
            velocityControl(i) = controlLimits.upper - 1e-03;
        }
    }

    return velocityControl;
}

RobotLibrary::Model::Limits
SerialLinkVelocityBase::compute_control_limits(const unsigned int &jointNumber)
{
    RobotLibrary::Model::Limits limits;

    double delta = _model->joint_positions()[jointNumber]
                 - _model->link(jointNumber)->joint().position_limits().lower;

    limits.lower = std::max(
        -delta * _controlFrequency,
        std::max(
            -_model->link(jointNumber)->joint().speed_limit(),
            -2 * std::sqrt(_maxJointAcceleration * delta)));

    delta = _model->link(jointNumber)->joint().position_limits().upper
          - _model->joint_positions()[jointNumber];

    limits.upper = std::min(
        delta * _controlFrequency,
        std::min(
            _model->link(jointNumber)->joint().speed_limit(),
            2 * std::sqrt(_maxJointAcceleration * delta)));

    if (limits.lower > limits.upper)
    {
        throw std::runtime_error(
            "[ERROR] [SERIAL LINK VELOCITY] compute_control_limits(): "
            "Lower limit for the '"
            + _model->link(jointNumber)->joint().name()
            + "' joint is greater than upper limit ("
            + std::to_string(limits.lower) + " > "
            + std::to_string(limits.upper) + ").");
    }

    return limits;
}

} } // namespace RobotLibrary::Control
