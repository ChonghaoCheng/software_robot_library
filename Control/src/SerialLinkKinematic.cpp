/**
 * @file    SerialLinkKinematic.cpp
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

#include <Control/SerialLinkKinematic.h>

namespace RobotLibrary { namespace Control {

  ///////////////////////////////////////////////////////////////////////////////////////////////////
 //                                          Constructor                                          //
///////////////////////////////////////////////////////////////////////////////////////////////////
SerialLinkKinematic::SerialLinkKinematic(std::shared_ptr<RobotLibrary::Model::KinematicTree> model,
		                                 const std::string &endpointName,
		                                 const RobotLibrary::Control::SerialLinkParameters &parameters)
: SerialLinkVelocityBase(model, endpointName, parameters)
{
    std::cout << "[INFO] [SERIAL LINK KINEMATICS] ";
    std::cout << "Performing VELOCITY control on the " + _model->name() + " robot.";
}
		                       
  ///////////////////////////////////////////////////////////////////////////////////////////////////
 //               Compute the endpoint velocity needed to track a given trajectory                //
///////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::VectorXd
SerialLinkKinematic::track_endpoint_trajectory(const RobotLibrary::Model::Pose &desiredPose,
                                               const Eigen::Vector<double,6>   &desiredVelocity,
                                               const Eigen::Vector<double,6>   &desiredAcceleration)
{
    (void)desiredAcceleration;                                                                      // Not needed in velocity control

    // NOTE: This method saves the magnitude of position and orientation error internally,
    //       so it can be queried after for analysing performance
    Eigen::Vector<double,6> poseError = pose_error(desiredPose);
    
	return resolve_endpoint_motion(desiredVelocity                                                  // Feedforward term
	                             + _cartesianPoseGain * poseError);                                 // Feedback term
}

} } // namespace
