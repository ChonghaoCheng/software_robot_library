/**
 * @file    SerialLinkVelocityBase.h
 * @brief   Shared resolved-rate control for serial link velocity controllers.
 */

#ifndef SERIAL_LINK_VELOCITY_BASE_H
#define SERIAL_LINK_VELOCITY_BASE_H

#include <Control/Core/SerialLinkBase.h>

namespace RobotLibrary { namespace Control {

/**
 * @brief Common twist-to-joint-velocity layer for velocity-level controllers.
 *
 * Derived controllers only need to decide which endpoint twist should be
 * commanded. This class performs the shared resolved-rate inverse kinematics,
 * joint-limit handling, singularity avoidance, redundancy handling, and joint
 * trajectory feedback using the state cached by the same controller object.
 */
class SerialLinkVelocityBase : public SerialLinkBase
{
    public:
        SerialLinkVelocityBase(
            std::shared_ptr<RobotLibrary::Model::KinematicTree> model,
            const std::string &endpointName,
            const RobotLibrary::Control::SerialLinkParameters &parameters = SerialLinkParameters());

        Eigen::VectorXd
        resolve_endpoint_motion(const Eigen::Vector<double,6> &endpointMotion) override;

        Eigen::VectorXd
        resolve_endpoint_twist(const Eigen::Vector<double,6> &twist) override;

        Eigen::VectorXd
        track_joint_trajectory(const Eigen::VectorXd &desiredPosition,
                               const Eigen::VectorXd &desiredVelocity,
                               const Eigen::VectorXd &desiredAcceleration) override;

    protected:
        RobotLibrary::Model::Limits
        compute_control_limits(const unsigned int &jointNumber) override;
};

} } // namespace RobotLibrary::Control

#endif
