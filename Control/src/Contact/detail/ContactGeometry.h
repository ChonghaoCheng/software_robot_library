/**
 * @file ContactGeometry.h
 * @brief Pure geometry helpers shared by moving-frame contact control and tests.
 */

#ifndef CONTROL_CONTACT_DETAIL_CONTACT_GEOMETRY_H
#define CONTROL_CONTACT_DETAIL_CONTACT_GEOMETRY_H

#include <Math/MathFunctions.h>
#include <Model/Pose.h>

#include <Eigen/Dense>

#include <algorithm>
#include <stdexcept>

namespace RobotLibrary { namespace Control { namespace detail {

inline Eigen::Vector3d
inward_contact_normal(const Eigen::Matrix3d &boardRotation,
                      const Eigen::Vector3d &normalAxisInBoard)
{
    const Eigen::Vector3d normal = boardRotation * normalAxisInBoard;
    if(!normal.allFinite() || normal.norm() <= 1e-9)
    {
        throw std::invalid_argument("Contact normal must be finite and nonzero.");
    }
    return normal.normalized();
}

inline double
signed_normal_coordinate(const Eigen::Vector3d &inwardNormal,
                         const Eigen::Vector3d &toolPosition,
                         const Eigen::Vector3d &boardPosition)
{
    return inwardNormal.dot(toolPosition - boardPosition);
}

struct AffineNormalForce
{
    Eigen::RowVectorXd map;
    double constant = 0.0;
};

/**
 * F_k = constant + map * U for p_tool,k = p_tool,0 + positionMap * U.
 */
inline AffineNormalForce
affine_normal_force(const double measuredNormalForce,
                    const double forceResponseGain,
                    const double currentSignedNormalCoordinate,
                    const Eigen::Vector3d &stageInwardNormal,
                    const Eigen::Vector3d &currentToolPosition,
                    const Eigen::Vector3d &stageBoardPosition,
                    const Eigen::Ref<const Eigen::MatrixXd> &positionMap)
{
    AffineNormalForce force;
    force.map = forceResponseGain * stageInwardNormal.transpose() * positionMap;
    force.constant = measuredNormalForce + forceResponseGain *
        (stageInwardNormal.dot(currentToolPosition - stageBoardPosition)
         - currentSignedNormalCoordinate);
    return force;
}

inline double
normal_force_from_coordinates(const double measuredNormalForce,
                              const double forceResponseGain,
                              const double currentSignedNormalCoordinate,
                              const double stageSignedNormalCoordinate)
{
    return measuredNormalForce
         + forceResponseGain * (stageSignedNormalCoordinate - currentSignedNormalCoordinate);
}

inline Eigen::Vector3d
local_orientation_reference(const Eigen::Quaterniond &currentOrientation,
                            const Eigen::Vector3d &currentRotationVector,
                            const Eigen::Quaterniond &desiredOrientation)
{
    return currentRotationVector
         + RobotLibrary::Math::quaternion_orientation_error(currentOrientation,
                                                            desiredOrientation);
}

inline Eigen::Vector<double,6>
impedance_pose_error(const RobotLibrary::Model::Pose &currentPose,
                     const RobotLibrary::Model::Pose &desiredPose)
{
    Eigen::Vector<double,6> error;
    error.head<3>() = desiredPose.translation() - currentPose.translation();
    error.tail<3>() = RobotLibrary::Math::quaternion_orientation_error(
        currentPose.quaternion(), desiredPose.quaternion());
    return error;
}

inline RobotLibrary::Model::Pose
propagate_rigid_frame(const RobotLibrary::Model::Pose &samplePose,
                      const Eigen::Vector3d &linearVelocity,
                      const Eigen::Vector3d &linearAcceleration,
                      const Eigen::Vector3d &spatialAngularVelocity,
                      const double elapsedTime)
{
    const double tau = std::max(0.0, elapsedTime);
    const Eigen::Vector3d position = samplePose.translation()
        + linearVelocity * tau + 0.5 * linearAcceleration * tau * tau;
    const Eigen::Matrix3d rotation =
        RobotLibrary::Math::so3_exponential(spatialAngularVelocity * tau)
        * samplePose.rotation();
    return RobotLibrary::Model::Pose(position, Eigen::Quaterniond(rotation).normalized());
}

} } } // namespace RobotLibrary::Control::detail

#endif
