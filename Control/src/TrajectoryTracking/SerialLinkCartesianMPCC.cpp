/**
 * @file    SerialLinkCartesianMPCC.cpp
 * @brief   Cartesian position MPCC with independent rotation-error tracking.
 */

#include <Control/TrajectoryTracking/SerialLinkCartesianMPCC.h>
#include <Math/MathFunctions.h>

#include <algorithm>

namespace RobotLibrary { namespace Control {

Eigen::Vector<double,6>
SerialLinkCartesianMPCC::cartesian_path_tangent(
    RobotLibrary::Trajectory::CartesianSpline &trajectory,
    const double progress)
{
    constexpr double step = 1e-3;
    const double s0 = std::clamp(progress, 0.0, 1.0);
    const double sForward = std::clamp(s0 + step, 0.0, 1.0);
    const double sPrevious = std::clamp(s0 - step, 0.0, 1.0);

    double denominator = sForward - s0;
    const RobotLibrary::Model::Pose referencePose =
        trajectory.pose_at_progress(s0);
    RobotLibrary::Model::Pose pose0 = referencePose;
    RobotLibrary::Model::Pose pose1;

    if(denominator > 1e-9)
    {
        pose1 = trajectory.pose_at_progress(sForward);
    }
    else
    {
        denominator = s0 - sPrevious;
        pose1 = pose0;
        pose0 = trajectory.pose_at_progress(sPrevious);
    }

    Eigen::Vector<double,6> tangent = Eigen::Vector<double,6>::Zero();
    if(denominator > 1e-9)
    {
        const Eigen::Vector3d positionTangent =
            (pose1.translation() - pose0.translation()) / denominator;
        // Express the parent-frame Cartesian tangent in the active path pose.
        // The same expression is obtained after any rigid left multiplication
        // by a board transform, so it remains compatible with SerialLinkMPCC.
        tangent.head<3>() =
            referencePose.quaternion().toRotationMatrix().transpose() * positionTangent;
    }

    // Rotation is not part of the Cartesian contour/lag projection, but its
    // path derivative is mandatory in e_R_dot = omega - tau_R(s) * sdot.
    // Without it, a continuously rotating reference can only be followed after
    // accumulating a feedback error of approximately |omega_ref| / K_R.
    tangent.tail<3>() = RobotLibrary::Math::so3_logarithm(
        pose0.quaternion().toRotationMatrix().transpose()
        * pose1.quaternion().toRotationMatrix()) / denominator;

    return tangent;
}

Eigen::Vector<double,6>
SerialLinkCartesianMPCC::path_tangent_at_progress(
    const double progress,
    const Eigen::Matrix3d &referenceRotation)
{
    const Eigen::Matrix3d frameRotation =
        trajectory_frame().transformInBase.block<3,3>(0,0);
    const Eigen::Matrix3d stageRotation =
        frameRotation
        * reference_trajectory().pose_at_progress(progress)
              .quaternion().toRotationMatrix();
    return mpcc_express_body_tangent_in_prediction_frame(
        cartesian_path_tangent(reference_trajectory(), progress),
        stageRotation,
        referenceRotation);
}

} } // namespace RobotLibrary::Control
