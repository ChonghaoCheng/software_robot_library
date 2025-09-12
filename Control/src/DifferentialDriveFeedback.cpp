/**
 * @file    DifferentialDriveFeedback.cpp
 * @author  Jon Woolfrey
 * @email   jonathan.woolfrey@gmail.com
 * @date    May 2025
 * @version 1.0
 * @brief   Source files for the DifferentialDriveFeedback control class.
 *
 * @details This class provides method for implementing nonlinear feedback control of a differential
 *          drive robot.
 * 
 * @copyright (c) 2025 Jon Woolfrey
 *
 * @license   OSCL - Free for non-commercial open-source use only.
 *            Commercial use requires a license.
 * 
 * @see https://github.com/Woolfrey/software_robot_library for more information.
 */
 
#include <Control/DifferentialDriveFeedback.h>

namespace RobotLibrary { namespace Control {

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                        Constructor                                             //
////////////////////////////////////////////////////////////////////////////////////////////////////
DifferentialDriveFeedback::DifferentialDriveFeedback(const RobotLibrary::Model::DifferentialDriveParameters &modelParameters,
                                                     const RobotLibrary::Control::DifferentialDriveFeedbackParameters &controlParameters)
: DifferentialDriveBase(controlParameters.controlFrequency,
                        controlParameters.minimumSafeDistance,
                        modelParameters,
                        controlParameters.qpsolver),
  _orientationGain(controlParameters.orientationGain),
  _xPositionGain(controlParameters.xPositionGain),
  _yPositionGain(controlParameters.yPositionGain)
{
    if (_xPositionGain   <= 0
    or  _yPositionGain   <= 0
    or  _orientationGain <= 0)
    {
        throw std::invalid_argument("[ERROR] [DIFFERENTIAL DRIVE FEEDBACK] Constructor: "
                                    "Feedback control gains must be positive, but "
                                    "the x position gain was " + std::to_string(_xPositionGain) + ", "
                                    "the y position gain was " + std::to_string(_yPositionGain) + ", and "
                                    "the orientation gain was " + std::to_string(_orientationGain) + ".");
    }
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                         Solve the control to track a desired trajectory                        //
////////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Vector2d
DifferentialDriveFeedback::track_trajectory(const RobotLibrary::Model::Pose2D &desiredPose,
                                            const Eigen::Vector2d &desiredVelocity,
                                            const std::vector<RobotLibrary::Model::Obstacle2D> &obstacles)
{
    // Kanayama, Y., Kimura, Y., Miyazaki, F., & Noguchi, T. (1990, May).
    // A stable tracking control method for an autonomous mobile robot.
    // In Proceedings., IEEE International Conference on Robotics and Automation (pp. 384-389). IEEE.
    
    // Pose error in the GLOBAL frame
    Eigen::Vector3d e = _pose.error(desiredPose);
    
    // Position error in the LOCAL frame
    double epsilon_x =  e[0] * cos(_pose.angle()) + e[1] * sin(_pose.angle());
    double epsilon_y = -e[0] * sin(_pose.angle()) + e[1] * cos(_pose.angle());
    
    /* Solve a QP problem of the form:
     *  min_u 1/2 (u_d - u)^T M (u_d - u)
     *  subject to: B * u <= z
     * Hessian H == M, and f == - M * u_d
     */
    
    Eigen::Matrix2d H;
    H << _mass,      0.0,
           0.0, _inertia;
                           
    Eigen::Vector2d f = { -_mass *    (desiredVelocity[0] * cos(e[2]) + _xPositionGain * epsilon_x), 
                          -_inertia * (desiredVelocity[1] + desiredVelocity[0] * epsilon_y + _yPositionGain * epsilon_y + _orientationGain * sin(e[2])) };
    
    RobotLibrary::Model::Limits linear, angular;
    
    compute_control_limits(linear, angular, velocity());

    _controlConstraintVector <<  linear.upper,
                                angular.upper,
                                -linear.lower,
                               -angular.lower;
     
    int n = obstacles.size();                                                                       // Makes referencing easier                          
                               
    _constraintMatrix.resize(4+n,2);
    _constraintMatrix.block(0,0,4,2) = _controlConstraintMatrix;
    _constraintMatrix.block(4,0,n,2) = _obstacleConstraintMatrix;
    
    _constraintVector.resize(4+n);
    _constraintVector.head(4) = _controlConstraintVector;
    _constraintVector.tail(n) = _obstacleConstraintVector;
    
    // NOTE: For some reason, in low dimensions, the QP solver needs a tiny kick-start
    Eigen::Vector2d u = solve(H,f,_constraintMatrix, _constraintVector, velocity());
    
    return {u[0], u[1]};
}

} } // Namespace                                                                                      
