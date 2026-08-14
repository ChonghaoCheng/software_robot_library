/**
 * @file    CartesianSpline.cpp
 * @author  Jon Woolfrey
 * @email   jonathan.woolfrey@gmail.com
 * @date    July 2025
 * @version 1.1
 * @brief   A class that defines trajectories for position & orientation in 3D space.
 * 
 * @details This class generates trajectories for position & orientation by using a spline to
 *          interpolate over a series of poses.
 * 
 * @copyright (c) 2025 Jon Woolfrey
 *
 * @license   OSCL - Free for non-commercial open-source use only.
 *            Commercial use requires a license.
 * 
 * @see https://github.com/Woolfrey/software_robot_library for more information.
 */
 
#include <Trajectory/CartesianSpline.h>
#include <Math/MathFunctions.h>

#include <algorithm>

namespace RobotLibrary { namespace Trajectory {

  //////////////////////////////////////////////////////////////////////////////////////////////////// 
 //                              Constructor for cubic splines                                     //
////////////////////////////////////////////////////////////////////////////////////////////////////
CartesianSpline::CartesianSpline(const std::vector<RobotLibrary::Model::Pose> &poses,
                                 const std::vector<double> &times,
                                 const Eigen::Vector<double,6> &startTwist)
{
    if(poses.size() < 2)
    {
        throw std::invalid_argument("[ERROR] [CARTESIAN SPLINE] Constructor: "
                                    "A minimum number of 2 poses is required to create a trajectory.");
    }
    else if(poses.size() != times.size())
    {
        throw std::invalid_argument("[ERROR] [CARTESIAN SPLINE] Constructor: "
                                    "Dimensions of arguments do not match. "
                                    "There were " + std::to_string(poses.size()) + " poses, but "
                                    + std::to_string(times.size()) + " times.");
    }
    
    // Convert the orientation from a quaternion to a vector
    std::vector<Eigen::VectorXd> positions(poses.size());                                           // NOTE: The SplineTrajectory class expects a Dynamic size vector
    for(int i = 0; i < positions.size(); i++)
    {   
        positions[i].resize(6);                                                                     // 3 for position, 3 for orientation
        
        positions[i].head(3) = poses[i].translation();                                              // First part is just the translation
        
        double angle = 2.0 * acos(std::clamp(poses[i].quaternion().normalized().w(), -1.0, 1.0));   // Get angle component
        
        if (abs(angle) < 1e-03) positions[i].tail(3).setZero();                                     // Practically zero
        else                    positions[i].tail(3) = angle * poses[i].quaternion().vec().normalized();
    }
    
    Eigen::Vector<double,6> startCoordinateVelocity = startTwist;
    startCoordinateVelocity.tail<3>() =
        RobotLibrary::Math::so3_left_jacobian_inverse(positions.front().tail<3>())
        * startTwist.tail<3>();

    this->_spline = SplineTrajectory(positions,times,startCoordinateVelocity);                      // Internal angular coordinates are rotation-vector derivatives.
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                            Get the desired state for the given time                            //
////////////////////////////////////////////////////////////////////////////////////////////////////
RobotLibrary::Trajectory::CartesianState
CartesianSpline::query_state(const double &time)
{
    RobotLibrary::Trajectory::State state = this->_spline.query_state(time);                        // Get the state as a 6x1 vector over real numbers
    
    // Now we need to convert the "position" vector to a pose

    double angle = state.position.tail(3).norm();                                                   // Norm of the angle * saxis component

    RobotLibrary::Model::Pose pose;                                                                 // We need to compute this

    if (angle < 1e-04) pose = RobotLibrary::Model::Pose(state.position.head(3),
                                                        Eigen::Quaterniond(1,0,0,0));               // Assume zero rotation
    else
    {
      Eigen::Vector<double,3> axis = state.position.tail(3).normalized();                           // Ensure magnitude of 1 
      
      pose = RobotLibrary::Model::Pose(state.position.head(3), 
                                       Eigen::Quaterniond(cos(0.5*angle),
                                                          sin(0.5*angle)*axis(0),
                                                          sin(0.5*angle)*axis(1),
                                                          sin(0.5*angle)*axis(2)));
    }

    const Eigen::Vector3d rotationVector = state.position.tail<3>();
    const Eigen::Vector3d rotationVectorRate = state.velocity.tail<3>();
    const Eigen::Matrix3d leftJacobian =
        RobotLibrary::Math::so3_left_jacobian(rotationVector);

    Eigen::Vector<double,6> physicalTwist = state.velocity;
    physicalTwist.tail<3>() = leftJacobian * rotationVectorRate;

    // alpha = J_l(r) rddot + Jdot_l(r, rdot) rdot. A centred directional
    // derivative keeps this consistent with the existing SO(3) Jacobian implementation.
    constexpr double derivativeStep = 1e-6;
    const Eigen::Matrix3d jacobianPlus = RobotLibrary::Math::so3_left_jacobian(
        rotationVector + derivativeStep * rotationVectorRate);
    const Eigen::Matrix3d jacobianMinus = RobotLibrary::Math::so3_left_jacobian(
        rotationVector - derivativeStep * rotationVectorRate);
    const Eigen::Matrix3d leftJacobianRate =
        (jacobianPlus - jacobianMinus) / (2.0 * derivativeStep);

    Eigen::Vector<double,6> physicalAcceleration = state.acceleration;
    physicalAcceleration.tail<3>() =
        leftJacobian * state.acceleration.tail<3>()
        + leftJacobianRate * rotationVectorRate;

    RobotLibrary::Trajectory::CartesianState returnValue =
        {pose, physicalTwist, physicalAcceleration};

    return returnValue;
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                       Map normalised progress s in [0,1] to trajectory time                     //
////////////////////////////////////////////////////////////////////////////////////////////////////
double
CartesianSpline::progress_to_time(const double &progress) const
{
    const double s = std::clamp(progress, 0.0, 1.0);
    return start_time() + s * (end_time() - start_time());                                           // Linear mapping
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                             Reference pose at a normalised progress                             //
////////////////////////////////////////////////////////////////////////////////////////////////////
RobotLibrary::Model::Pose
CartesianSpline::pose_at_progress(const double &progress)
{
    return this->query_state(progress_to_time(progress)).pose;
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                       Body-frame SE(3) tangent w.r.t. normalised progress                       //
////////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Matrix<double,6,1>
CartesianSpline::tangent_at_progress(const double &progress, const double &step)
{
    const double s0    = std::clamp(progress, 0.0, 1.0);
    const double s1    = std::clamp(s0 + step, 0.0, 1.0);
    const double sPrev = std::clamp(s0 - step, 0.0, 1.0);

    double denominator = s1 - s0;
    Eigen::Matrix4d T0 = pose_at_progress(s0).as_matrix();
    Eigen::Matrix4d T1;
    if(denominator > 1e-9)
    {
        T1 = pose_at_progress(s1).as_matrix();                                                       // Forward difference
    }
    else
    {
        denominator = s0 - sPrev;                                                                    // At the upper end, use a backward difference
        T1 = T0;
        T0 = pose_at_progress(sPrev).as_matrix();
    }

    if(denominator <= 1e-9) return Eigen::Matrix<double,6,1>::Zero();

    return RobotLibrary::Math::se3_logarithm(RobotLibrary::Math::se3_inverse(T0) * T1) / denominator;
}

} } // namespace
