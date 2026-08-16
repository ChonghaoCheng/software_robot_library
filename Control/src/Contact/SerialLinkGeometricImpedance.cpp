/**
 * @file    SerialLinkGeometricImpedance.cpp
 * @author  ChonghaoCheng
 * @date    June 2026
 * @version 1.0.0
 *
 * @brief   Geometric (SE(3)) impedance control with force feedforward / tracking.
 *
 * @see     SerialLinkGeometricImpedance.h
 * @see     Seo et al., arXiv:2211.07945
 *
 * @license OSCL - Free for non-commercial open-source use only.
 */

#include <Control/Contact/SerialLinkGeometricImpedance.h>

namespace RobotLibrary { namespace Control {

  ///////////////////////////////////////////////////////////////////////////////////////////////////
 //                                            Constructor                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////
SerialLinkGeometricImpedance::SerialLinkGeometricImpedance(std::shared_ptr<RobotLibrary::Model::KinematicTree> model,
                                                           const std::string &endpointName,
                                                           const RobotLibrary::Control::SerialLinkParameters &parameters)
: SerialLinkBase(model, endpointName, parameters)
{
    unsigned int n = _model->number_of_joints();

    _desiredConfiguration.resize(n);

    for (int i = 0; i < n; ++i)
    {
        RobotLibrary::Model::Limits limits = _model->link(i)->joint().position_limits();

        _desiredConfiguration(i) = (limits.lower + limits.upper) / 2.0;
    }

    std::cout << "[INFO] [SERIAL LINK GEOMETRIC IMPEDANCE] "
              << "Performing TORQUE control on the " + _model->name() + " robot.\n";

    std::cout << "[INFO] [SERIAL LINK GEOMETRIC IMPEDANCE] "
              << "Gravity compensation is " << (_gravityCompensation ? "ENABLED" : "DISABLED")
              << "; force-tracking gains default to zero (pure feedforward)." << std::endl;
}

  ///////////////////////////////////////////////////////////////////////////////////////////////////
 //                            Track a Cartesian trajectory under contact                         //
///////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::VectorXd
SerialLinkGeometricImpedance::track_endpoint_trajectory(const RobotLibrary::Model::Pose &desiredPose,
                                                        const Eigen::Vector<double,6>   &desiredVelocity,
                                                        const Eigen::Vector<double,6>   &desiredAcceleration)
{
    // Geometrically consistent SE(3) pose error (position + angle*axis orientation, base frame).
    // NOTE: pose_error() also stores the position/orientation error magnitudes internally.
    Eigen::Vector<double,6> poseError = pose_error(desiredPose);

    Eigen::Vector<double,6> impedanceWrench = _cartesianPoseGain     * poseError
                                            + _cartesianVelocityGain * (desiredVelocity - endpoint_velocity());

    return map_wrench_to_torque(impedanceWrench + force_wrench());
}

  ///////////////////////////////////////////////////////////////////////////////////////////////////
 //                       Damping for a desired twist plus the feedforward wrench                  //
///////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::VectorXd
SerialLinkGeometricImpedance::resolve_endpoint_twist(const Eigen::Vector<double,6> &twist)
{
    Eigen::Vector<double,6> impedanceWrench = _cartesianVelocityGain * (twist - endpoint_velocity());

    return map_wrench_to_torque(impedanceWrench + force_wrench());
}

  ///////////////////////////////////////////////////////////////////////////////////////////////////
 //                            Compute the force feedforward + tracking term                      //
///////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Vector<double,6>
SerialLinkGeometricImpedance::force_wrench()
{
    Eigen::Vector<double,6> forceError = _desiredWrench - _measuredWrench;

    // Integrate the wrench error (with simple anti-windup) only if integral action is used.
    if (!_forceIntegralGain.isZero())
    {
        _forceIntegral += forceError / _controlFrequency;

        for (int i = 0; i < 6; ++i)
        {
            _forceIntegral(i) = std::clamp(_forceIntegral(i), -_forceIntegralLimit, _forceIntegralLimit);
        }
    }

    return _desiredWrench + _forceGain * forceError + _forceIntegralGain * _forceIntegral;
}

  ///////////////////////////////////////////////////////////////////////////////////////////////////
 //                Map an end-effector wrench to joint torques (+ null space + gravity)           //
///////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::VectorXd
SerialLinkGeometricImpedance::map_wrench_to_torque(const Eigen::Vector<double,6> &wrench)
{
    using namespace Eigen;

    unsigned int n = _model->number_of_joints();

    VectorXd nullSpaceTask(n);

    for (int i = 0; i < n; ++i)
    {
        nullSpaceTask[i] = _jointPositionGains[i]/5.0 * (_desiredConfiguration[i] - _model->joint_positions()[i])
                         - _jointVelocityGains[i]/5.0 *  _model->joint_velocities()[i];
    }

    const auto &[Q, R] = RobotLibrary::Math::schwarz_rutishauser(_jacobianMatrix.transpose(), 1e-06);

    VectorXd jointTorques = _jacobianMatrix.transpose() * wrench
                          + (MatrixXd::Identity(n,n) - Q * Q.transpose()) * nullSpaceTask;

    if (_gravityCompensation) jointTorques += _model->joint_gravity_vector();                        // g(q) holds the arm against gravity

    // Ensure the joint torque limits are obeyed
    for (int i = 0; i < n; ++i)
    {
        const auto &limits = compute_control_limits(i);

        jointTorques[i] = std::clamp(jointTorques[i], limits.lower + 1e-03, limits.upper - 1e-03);
    }

    return jointTorques;
}

  ///////////////////////////////////////////////////////////////////////////////////////////////////
 //                                       Not used in this class                                  //
///////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::VectorXd
SerialLinkGeometricImpedance::resolve_endpoint_motion(const Eigen::Vector<double,6> &endpointMotion)
{
    throw std::runtime_error("[ERROR] [SERIAL LINK GEOMETRIC IMPEDANCE] resolve_endpoint_motion(): "
                             "This method is not used in this class. Did you mean to call `map_wrench_to_torque()'?");

    return Eigen::VectorXd::Zero(_model->number_of_joints());
}

  ///////////////////////////////////////////////////////////////////////////////////////////////////
 //                              Joint-space PD with gravity compensation                         //
///////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::VectorXd
SerialLinkGeometricImpedance::track_joint_trajectory(const Eigen::VectorXd &desiredPosition,
                                                     const Eigen::VectorXd &desiredVelocity,
                                                     const Eigen::VectorXd &desiredAcceleration)
{
    unsigned int n = _model->number_of_joints();

    Eigen::VectorXd jointTorques(n);

    Eigen::VectorXd gravity = _gravityCompensation ? _model->joint_gravity_vector()
                                                   : Eigen::VectorXd::Zero(n);

    for (int i = 0; i < n; ++i)
    {
        RobotLibrary::Model::Limits positionLimits = _model->link(i)->joint().position_limits();
        double speedLimit = _model->link(i)->joint().speed_limit();

        double positionReference = std::clamp(desiredPosition[i], positionLimits.lower + 1e-03, positionLimits.upper - 1e-03);

        double velocityReference = std::clamp(desiredVelocity[i], -speedLimit + 1e-03, speedLimit - 1e-03);

        jointTorques[i] = _jointVelocityGains[i] * (velocityReference - _model->joint_velocities()[i])
                        + _jointPositionGains[i] * (positionReference - _model->joint_positions()[i])
                        + gravity[i];

        RobotLibrary::Model::Limits torqueLimits = compute_control_limits(i);

        jointTorques[i] = std::clamp(jointTorques[i], torqueLimits.lower + 1e-03, torqueLimits.upper - 1e-03);
    }

    return jointTorques;
}

  ///////////////////////////////////////////////////////////////////////////////////////////////////
 //                                  Set the force-tracking gains                                 //
///////////////////////////////////////////////////////////////////////////////////////////////////
void
SerialLinkGeometricImpedance::set_force_gains(const Eigen::Matrix<double,6,6> &proportional,
                                             const Eigen::Matrix<double,6,6> &integral)
{
    _forceGain         = proportional;
    _forceIntegralGain = integral;
    _forceIntegral.setZero();
}

  ///////////////////////////////////////////////////////////////////////////////////////////////////
 //                      Set the desired joint configuration for redundancy                       //
///////////////////////////////////////////////////////////////////////////////////////////////////
bool
SerialLinkGeometricImpedance::set_desired_configuration(const Eigen::VectorXd &configuration)
{
    if (configuration.size() != _model->number_of_joints())
    {
        std::cout << "[ERROR] [SERIAL LINK GEOMETRIC IMPEDANCE] set_desired_configuration(): "
                  << "Input argument had " << configuration.size() << " elements "
                  << "but expected " << _model->number_of_joints() << "." << std::endl;

        return false;
    }

    for (int i = 0; i < _model->number_of_joints(); ++i)
    {
        RobotLibrary::Model::Limits limits = _model->link(i)->joint().position_limits();

        _desiredConfiguration(i) = std::clamp(configuration(i), limits.lower + 1e-03, limits.upper - 1e-03);
    }

    return true;
}

  ///////////////////////////////////////////////////////////////////////////////////////////////////
 //                                  Compute the control limits                                   //
///////////////////////////////////////////////////////////////////////////////////////////////////
RobotLibrary::Model::Limits
SerialLinkGeometricImpedance::compute_control_limits(const unsigned int &jointNumber)
{
    RobotLibrary::Model::Limits limits;

    limits.lower = -_model->link(jointNumber)->joint().effort_limit();
    limits.upper =  _model->link(jointNumber)->joint().effort_limit();

    if (limits.lower > limits.upper)
    {
        throw std::logic_error("[ERROR] [SERIAL LINK GEOMETRIC IMPEDANCE] compute_control_limits(): "
                               "Lower limit for the '" + _model->link(jointNumber)->joint().name() + "' joint is greater than "
                               "upper limit (" + std::to_string(limits.lower) + " > " + std::to_string(limits.upper) + ").");
    }

    return limits;
}

} } // namespace
