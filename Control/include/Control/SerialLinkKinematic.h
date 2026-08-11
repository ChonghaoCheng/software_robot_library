/**
 * @file    SerialLinkKinematic.h
 * @author  Jon Woolfrey
 * @email   jonathan.woolfrey@gmail.com
 * @date    August 2025
 * @version 2.0.1
 *
 * @brief   Computes velocity (position) feedback control for a serial link robot arm.
 * 
 * @details This class contains methods for performing velocity control of a serial link robot arm
 *          in both Cartesian and joint space. The fundamental feedforward + feedback control is given by:
 *          control velocity = desired velocity + gain * (desired position - actual position).
 * 
 * @copyright (c) 2025 Jon Woolfrey
 *
 * @license   OSCL - Free for non-commercial open-source use only.
 *            Commercial use requires a license.
 *
 * @see https://github.com/Woolfrey/software_robot_library for more information.
 * @see https://github.com/Woolfrey/software_simple_qp for the optimisation algorithm used in the control.
 */

#ifndef SERIAL_LINK_KINEMATIC_H
#define SERIAL_LINK_KINEMATIC_H

#include <Control/SerialLinkVelocityBase.h>

#include <memory>

namespace RobotLibrary { namespace Control {

/**
 * @brief Algorithms for velocity control of a serial link robot arm.
 */
class SerialLinkKinematic : public SerialLinkVelocityBase
{
	public:
		/**
		 * @brief Constructor.
		 * @param model A pointer to a KinematicTree object.
		 * @param endpointName The name of the reference frame in the KinematicTree to be controlled.
		 */
		SerialLinkKinematic(std::shared_ptr<RobotLibrary::Model::KinematicTree> model,
		                    const std::string &endpointName,
		                    const RobotLibrary::Control::SerialLinkParameters &parameters = SerialLinkParameters());
		
		/**
		 * @brief Solve the joint velocities required to track a Cartesian trajectory.
		 * @param desiredPose The desired position & orientation (pose) for the endpoint.
		 * @param desiredVel The desired linear & angular velocity (twist) for the endpoint.
		 * @param desiredAcc Not used in velocity control.
		 * @return The joint velocities (nx1) required to track the trajectory.
		 */
		Eigen::VectorXd
		track_endpoint_trajectory(const RobotLibrary::Model::Pose &desiredPose,
					              const Eigen::Vector<double,6>   &desiredVelocity,
					              const Eigen::Vector<double,6>   &desiredAcceleration) override;
	
};                                                                                                  // Semicolon needed after a class declaration

} } // namespace

#endif
