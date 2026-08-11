/**
 * @file    SerialLinkCartesianMPCC.cpp
 * @brief   Cartesian position MPCC with independent rotation-error tracking.
 */

#include <Control/SerialLinkCartesianMPCC.h>

#include <algorithm>

namespace RobotLibrary { namespace Control {

Eigen::Vector<double,6>
SerialLinkCartesianMPCC::path_tangent_at_progress(
    const double progress,
    const Eigen::Matrix3d &referenceRotation)
{
    constexpr double step = 1e-3;
    const double s0 = std::clamp(progress, 0.0, 1.0);
    const double sForward = std::clamp(s0 + step, 0.0, 1.0);
    const double sPrevious = std::clamp(s0 - step, 0.0, 1.0);

    double denominator = sForward - s0;
    RobotLibrary::Model::Pose pose0 = reference_trajectory().pose_at_progress(s0);
    RobotLibrary::Model::Pose pose1;

    if(denominator > 1e-9)
    {
        pose1 = reference_trajectory().pose_at_progress(sForward);
    }
    else
    {
        denominator = s0 - sPrevious;
        pose1 = pose0;
        pose0 = reference_trajectory().pose_at_progress(sPrevious);
    }

    Eigen::Vector<double,6> tangent = Eigen::Vector<double,6>::Zero();
    if(denominator > 1e-9)
    {
        const Eigen::Vector3d positionTangent =
            (pose1.translation() - pose0.translation()) / denominator;
        tangent.head<3>() = referenceRotation.transpose() * positionTangent;
    }

    return tangent;
}

} } // namespace RobotLibrary::Control
