#include <Trajectory/CartesianSpline.h>
#include <Trajectory/CartesianTrajectoryFrame.h>
#include <Math/MathFunctions.h>

#include <Eigen/Geometry>
#include <cmath>
#include <iostream>

namespace {

bool near(const Eigen::VectorXd &a, const Eigen::VectorXd &b, double tolerance)
{
    return a.size() == b.size() && (a - b).norm() <= tolerance;
}

Eigen::Matrix4d pose_matrix(const RobotLibrary::Model::Pose &pose)
{
    Eigen::Matrix4d transform = Eigen::Matrix4d::Identity();
    transform.block<3,3>(0,0) = pose.quaternion().toRotationMatrix();
    transform.block<3,1>(0,3) = pose.translation();
    return transform;
}

} // namespace

int main()
{
    using RobotLibrary::Trajectory::CartesianTrajectoryFrameState;

    CartesianTrajectoryFrameState frame;
    frame.transformInBase.block<3,3>(0,0) =
        Eigen::AngleAxisd(M_PI_2, Eigen::Vector3d::UnitZ()).toRotationMatrix();
    frame.transformInBase.block<3,1>(0,3) = Eigen::Vector3d(1.0, 2.0, 3.0);
    frame.twistInBase << 0.1, 0.2, 0.3, 0.0, 0.0, 2.0;

    const RobotLibrary::Model::Pose localPose(
        Eigen::Vector3d(0.5, 0.0, 0.0), Eigen::Quaterniond::Identity());
    const RobotLibrary::Model::Pose basePose =
        RobotLibrary::Trajectory::express_pose_in_base(frame, localPose);
    if(!near(basePose.translation(), Eigen::Vector3d(1.0, 2.5, 3.0), 1e-12)) return 1;

    Eigen::Vector<double,6> localTwist;
    localTwist << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0;
    const Eigen::Vector<double,6> baseTwist =
        RobotLibrary::Trajectory::express_point_twist_in_base(frame, localPose, localTwist);
    Eigen::Vector<double,6> expectedTwist;
    expectedTwist << -0.9, 1.2, 0.3, -1.0, 0.0, 2.0;
    if(!near(baseTwist, expectedTwist, 1e-12)) return 2;

    const Eigen::Vector3d r0(0.4, -0.2, 0.3);
    const Eigen::Vector3d r1(0.7, 0.1, 0.5);
    std::vector<RobotLibrary::Model::Pose> poses{
        {Eigen::Vector3d::Zero(), Eigen::Quaterniond(RobotLibrary::Math::so3_exponential(r0))},
        {Eigen::Vector3d(0.1, 0.2, 0.3), Eigen::Quaterniond(RobotLibrary::Math::so3_exponential(r1))}
    };
    Eigen::Vector<double,6> startTwist = Eigen::Vector<double,6>::Zero();
    startTwist.tail<3>() = Eigen::Vector3d(0.2, -0.1, 0.3);
    RobotLibrary::Trajectory::CartesianSpline spline(poses, {0.0, 2.0}, startTwist);

    const double time = 0.7;
    const double h = 1e-5;
    const auto state = spline.query_state(time);
    const Eigen::Matrix3d before = spline.query_state(time - h).pose.quaternion().toRotationMatrix();
    const Eigen::Matrix3d after = spline.query_state(time + h).pose.quaternion().toRotationMatrix();
    const Eigen::Vector3d finiteDifferenceOmega =
        RobotLibrary::Math::so3_logarithm(after * before.transpose()) / (2.0 * h);
    if(!near(state.twist.tail<3>(), finiteDifferenceOmega, 2e-5)) return 3;

    const auto beforeState = spline.query_state(time - h);
    const auto afterState = spline.query_state(time + h);
    const Eigen::Vector3d finiteDifferenceAlpha =
        (afterState.twist.tail<3>() - beforeState.twist.tail<3>()) / (2.0 * h);
    if(!near(state.acceleration.tail<3>(), finiteDifferenceAlpha, 2e-4)) return 4;

    const Eigen::Matrix4d local0 = pose_matrix(poses.front());
    const Eigen::Matrix4d local1 = pose_matrix(poses.back());
    const Eigen::Matrix4d left = frame.transformInBase;
    const Eigen::Matrix4d relativeLocal = local0.inverse() * local1;
    const Eigen::Matrix4d relativeBoard = (left * local0).inverse() * (left * local1);
    if((relativeLocal - relativeBoard).norm() > 1e-12) return 5;

    return 0;
}
