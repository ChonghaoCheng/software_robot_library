/**
 * @file    DifferentialDriveBase.h
 * @author  Jon Woolfrey
 * @email   jonathan.woolfrey@gmail.com
 * @date    May 2025
 * @version 1.0
 * @brief   A base class to standardise all control classes of differential drive robots.
 * 
 * @copyright (c) 2025 Jon Woolfrey
 *
 * @license   OSCL - Free for non-commercial open-source use only.
 *            Commercial use requires a license.
 * 
 * @see https://github.com/Woolfrey/software_robot_library for more information.
 */
 
#ifndef DIFFERENTIAL_DRIVE_BASE_H
#define DIFFERENTIAL_DRIVE_BASE_H

#include <Math/QPSolver.h>
#include <Model/DifferentialDrive.h>

namespace RobotLibrary { namespace Control {

class DifferentialDriveBase : public QPSolver<double>,
                              public RobotLibrary::Model::DifferentialDrive
{
    public:
    
        DifferentialDriveBase(const double &controlFrequency,
                              const double &minimumSafeDistance,
                              const RobotLibrary::Model::DifferentialDriveParameters &modelParameters,
                              const SolverOptions<double> &solverOptions);
                              
        /**
         * @brief Get the predicted pose at the next time step given the current state.
         */
        RobotLibrary::Model::Pose2D
        predicted_pose()
        { 
            return RobotLibrary::Model::DifferentialDrive::predicted_pose(_pose, velocity(), _controlFrequency);
        }
    
    protected:

        double _controlFrequency;                                                                   ///< Rate at which control commands are sent to robot
        
        double _minimumSafeDistance;                                                                ///< Extra safety distance used in collision avoidance

        Eigen::Matrix<double, Eigen::Dynamic, 2> _constraintMatrix;                                 ///< For the QP solver
        
        Eigen::Matrix<double, 4, 2> _controlConstraintMatrix;                                       ///< Limits on the instantaneous velocity
        
        Eigen::Matrix<double, Eigen::Dynamic, 2> _obstacleConstraintMatrix;                         ///< Part of the control barrier function
        
        Eigen::Vector<double, Eigen::Dynamic> _constraintVector;                                    ///< Limits on instantaneous velocity
        
        Eigen::Vector<double, 4> _controlConstraintVector;                                          ///< For the QP solver
        
        Eigen::Vector<double, Eigen::Dynamic> _obstacleConstraintVector;                            ///< Part of the control barrier function
        
        /**
         * @brief Compute the instantaneous limits on the linear & angular velocity.
         * @param Storage for the linear velocity limits.
         * @param Storage for the angular velocity limits.
         */
        void
        compute_control_limits(RobotLibrary::Model::Limits &linear,
                               RobotLibrary::Model::Limits &angular,
                               const Eigen::Vector2d &currentVelocity);
};
 
} } // namespace

#endif
