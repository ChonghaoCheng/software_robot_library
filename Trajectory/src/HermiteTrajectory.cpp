/**
 * @file    HermiteTrajectory.cpp
 * @author  Jon Woolfrey
 * @email   jonathan.woolfrey@gmail.com
 * @date    March 2026
 * @version 1.0
 * @brief   Source files for the HermiteTrajectory class.
 * 
 * @details This class uses a Hermite spline to generate a trajectory for a mobile robot whose
 *          kinematics is described by a unicycle model. Its configuration is SE(2), and its control
 *          space is linear velocity (m/s) and angular velocity (rad/s).
 *
 * @copyright (c) 2026 Jon Woolfrey
 *
 * @license   OSCL - Free for non-commercial open-source use only.
 *            Commercial use requires a license.
 * 
 * @see https://github.com/Woolfrey/software_robot_library for more information.
 */
 
#include <Trajectory/HermiteTrajectory.h>

namespace RobotLibrary { namespace Trajectory {

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                         Constructor                                            //
////////////////////////////////////////////////////////////////////////////////////////////////////
HermiteTrajectory::HermiteTrajectory(const RobotLibrary::Model::Pose2D &startPose,
                                     const RobotLibrary::Model::Pose2D &endPose,
                                     const double &startTime,
                                     const double &endTime,
                                     const Eigen::Vector3d &startTwist)
: _startPose(startPose),
  _endPose(endPose),
  _startTime(startTime),
  _endTime(endTime),
  _startTwist(startTwist)
{
    // Ensure inputs are logical
    if (startTime >= endTime)
    {
        throw std::invalid_argument("[ERROR] [UNICYCLE TRAJECTORY] Constructor: "
                                    "Start time " + std::to_string(startTime) + " >= "
                                    "end time " + std::to_string (endTime) + ".");
    }
    
    using namespace RobotLibrary::Math;                                                             // For brevity
    using namespace RobotLibrary::Trajectory;
    
    double angle     = startPose.angle();
    double velocity  = startTwist[0] * cos(angle) + startTwist[1] * sin(angle);
    double curvature = startTwist[2] / (velocity + 1e-12);
    
    // Set up the x-position spline
    FunctionPoint xStart;
    xStart.value            = startPose.translation()[0];
    xStart.firstDerivative  = cos(angle);
    xStart.secondDerivative = - curvature * sin(angle);
    
    FunctionPoint xEnd;
    xEnd.value            = endPose.translation()[0];
    xEnd.firstDerivative  = cos(endPose.angle());
    xEnd.secondDerivative = 0.0;
    
    _xSpline = Polynomial(xStart, xEnd, 0.0, 1.0, 5);                                               // Cubic polynomial
    
    // Set up the y-position spline
    FunctionPoint yStart;
    yStart.value            = startPose.translation()[1];
    yStart.firstDerivative  = sin(angle);
    yStart.secondDerivative = curvature * cos(angle);
    
    FunctionPoint yEnd;
    yEnd.value            = endPose.translation()[1];
    yEnd.firstDerivative  = sin(endPose.angle());
    yEnd.secondDerivative = 0.0;
    
    _ySpline = Polynomial(yStart, yEnd, 0.0, 1.0, 5);                                               // Cubic polynomial

    // Set up the arc length trajectory
    FunctionPoint arcParameterStart;
    arcParameterStart.value            = 0.0;
    arcParameterStart.firstDerivative  = velocity;
    arcParameterStart.secondDerivative = 0.0;
    
    FunctionPoint arcParameterEnd;
    arcParameterEnd.value            = 1.0;
    arcParameterEnd.firstDerivative  = 0.0;
    arcParameterEnd.secondDerivative = 0.0;
    
    _arcParameterTrajectory = Polynomial(arcParameterStart, arcParameterEnd, startTime, endTime, 5); // Quintic polynomial function of time?
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                               Get the desired state for the given time                         //
////////////////////////////////////////////////////////////////////////////////////////////////////
RobotLibrary::Trajectory::PlanarState
HermiteTrajectory::query_state(const double &time)
{
    using namespace RobotLibrary::Model;
    using namespace RobotLibrary::Trajectory;
    using namespace Eigen;
    
    PlanarState desiredState;                                                                       // Value to be returned
    
    if (time <= _startTime)
    {
        desiredState.pose         = _startPose;
        desiredState.twist        = _startTwist;
        desiredState.acceleration = {0.0, 0.0, 0.0};
    }
    else if (time > _endTime)
    {

        desiredState.pose         = _endPose;
        desiredState.twist        = {0.0, 0.0, 0.0};                                                // Assume zero velocity at the end
        desiredState.acceleration = {0.0, 0.0, 0.0};
    }
    else
    {
        const auto &[s, sdot, sddot] = _arcParameterTrajectory.evaluate_point(time);                // Arc length as a function of time
        const auto &[x, dx, ddx]     = _xSpline.evaluate_point(s);                                  // x-coordinate as a function of arc length
        const auto &[y, dy, ddy]     = _ySpline.evaluate_point(s);                                  // y-coordinate as a function of arc length

        double angle = std::atan2(dy,dx);                                                           // i.e heading angle
        
        double denominator = dx * dx + dy * dy + 1e-12;
        
        desiredState.pose = Pose2D(x, y, angle);
        
        desiredState.twist= {dx * sdot,
                             dy * sdot,
                             sdot * (dx * ddy - dy * ddx) / pow(dx * dx + dy * dy + 1e-12, 1.5)};
        
        desiredState.acceleration = {dx * sddot + ddx * sdot,
                                     dy * sddot + ddy * sdot,
                                     0.0}; // NOTE: This requires d^3 x/ds^3 which we don't currently have...
    }

    return desiredState;
}      

} } // namespace
