/**
 * @file    CartesianTrajectoryFrame.h
 * @brief   Parent-frame state and coordinate transforms for Cartesian trajectories.
 */

#ifndef CARTESIAN_TRAJECTORY_FRAME_H
#define CARTESIAN_TRAJECTORY_FRAME_H

#include <Model/Pose.h>
#include <Trajectory/DataStructures.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace RobotLibrary { namespace Trajectory {

/**
 * @brief Pose and point-twist of a Cartesian trajectory's parent frame.
 *
 * transformInBase is ^B T_F. twistInBase is [^B v_F; ^B omega_F], where
 * v_F is the velocity of the frame origin, both components expressed in B.
 */
struct CartesianTrajectoryFrameState
{
    Eigen::Matrix4d transformInBase = Eigen::Matrix4d::Identity();
    Eigen::Vector<double,6> twistInBase = Eigen::Vector<double,6>::Zero();
    /** Source measurement timestamp; NaN when the producer has no clock stamp. */
    double measurementTimeSeconds = std::numeric_limits<double>::quiet_NaN();
};

inline void
validate_trajectory_frame(const CartesianTrajectoryFrameState &frame)
{
    if(not frame.transformInBase.allFinite() or not frame.twistInBase.allFinite())
    {
        throw std::invalid_argument(
            "[ERROR] [CARTESIAN TRAJECTORY FRAME] Frame pose and twist must be finite.");
    }

    const Eigen::Matrix3d rotation = frame.transformInBase.block<3,3>(0,0);
    if((rotation.transpose() * rotation - Eigen::Matrix3d::Identity()).norm() > 1e-6
       or std::abs(rotation.determinant() - 1.0) > 1e-6
       or (frame.transformInBase.row(3) - Eigen::RowVector4d(0.0, 0.0, 0.0, 1.0)).norm() > 1e-9)
    {
        throw std::invalid_argument(
            "[ERROR] [CARTESIAN TRAJECTORY FRAME] transformInBase must be an SE(3) transform.");
    }
}

inline RobotLibrary::Model::Pose
express_pose_in_base(const CartesianTrajectoryFrameState &frame,
                     const RobotLibrary::Model::Pose &poseInFrame)
{
    const Eigen::Matrix3d rotation = frame.transformInBase.block<3,3>(0,0);
    const Eigen::Vector3d translation = frame.transformInBase.block<3,1>(0,3);
    const Eigen::Vector3d positionInBase = translation + rotation * poseInFrame.translation();
    const Eigen::Quaterniond orientationInBase =
        Eigen::Quaterniond(rotation) * poseInFrame.quaternion();
    return RobotLibrary::Model::Pose(positionInBase, orientationInBase.normalized());
}

inline Eigen::Vector<double,6>
express_point_twist_in_base(const CartesianTrajectoryFrameState &frame,
                            const RobotLibrary::Model::Pose &poseInFrame,
                            const Eigen::Vector<double,6> &twistInFrame)
{
    const Eigen::Matrix3d rotation = frame.transformInBase.block<3,3>(0,0);
    const Eigen::Vector3d offsetInBase = rotation * poseInFrame.translation();

    Eigen::Vector<double,6> twistInBase;
    twistInBase.head<3>() = frame.twistInBase.head<3>()
                            + frame.twistInBase.tail<3>().cross(offsetInBase)
                            + rotation * twistInFrame.head<3>();
    twistInBase.tail<3>() = frame.twistInBase.tail<3>()
                            + rotation * twistInFrame.tail<3>();
    return twistInBase;
}

inline CartesianState
express_state_in_base(const CartesianTrajectoryFrameState &frame,
                      const CartesianState &stateInFrame)
{
    CartesianState stateInBase;
    stateInBase.pose = express_pose_in_base(frame, stateInFrame.pose);
    stateInBase.twist = express_point_twist_in_base(frame, stateInFrame.pose, stateInFrame.twist);

    // The frame contract intentionally contains pose/twist only. Rotate the
    // relative acceleration; moving-frame transport acceleration remains a
    // responsibility of dedicated moving-frame controllers.
    const Eigen::Matrix3d rotation = frame.transformInBase.block<3,3>(0,0);
    stateInBase.acceleration.head<3>() = rotation * stateInFrame.acceleration.head<3>();
    stateInBase.acceleration.tail<3>() = rotation * stateInFrame.acceleration.tail<3>();
    return stateInBase;
}

} } // namespace RobotLibrary::Trajectory

#endif
