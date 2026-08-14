#include <Control/SerialLinkMPCC.h>
#include <Trajectory/CartesianSpline.h>

#include <Eigen/Geometry>
#include <iostream>

int main()
{
    using RobotLibrary::Model::Pose;
    using RobotLibrary::Trajectory::CartesianSpline;

    const Pose start(
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Quaterniond(Eigen::AngleAxisd(0.2, Eigen::Vector3d::UnitX())));
    const Pose finish(
        Eigen::Vector3d(0.35, -0.12, 0.18),
        Eigen::Quaterniond(
            Eigen::AngleAxisd(1.1, Eigen::Vector3d(0.2, 0.9, -0.1).normalized())));
    CartesianSpline trajectory(
        {start, finish}, {0.0, 2.0}, Eigen::Vector<double,6>::Zero());

    constexpr double predictionFrameProgress = 0.2;
    constexpr double stageProgress = 0.72;
    const Eigen::Matrix3d predictionRotation =
        trajectory.pose_at_progress(predictionFrameProgress)
            .quaternion().toRotationMatrix();
    const Eigen::Matrix3d stageRotation =
        trajectory.pose_at_progress(stageProgress)
            .quaternion().toRotationMatrix();
    const Eigen::Vector<double,6> stageBodyTangent =
        trajectory.tangent_at_progress(stageProgress);

    Eigen::Vector<double,6> expectedInPredictionFrame;
    const Eigen::Matrix3d stageToPrediction =
        predictionRotation.transpose() * stageRotation;
    expectedInPredictionFrame.head<3>() =
        stageToPrediction * stageBodyTangent.head<3>();
    expectedInPredictionFrame.tail<3>() =
        stageToPrediction * stageBodyTangent.tail<3>();

    const Eigen::Vector<double,6> implemented =
        RobotLibrary::Control::mpcc_express_body_tangent_in_prediction_frame(
            stageBodyTangent, stageRotation, predictionRotation);
    const double error = (implemented - expectedInPredictionFrame).norm();
    if(error > 1e-9)
    {
        std::cerr << "MPCC stage tangent changes coordinates across the horizon: "
                  << error << '\n';
        return 1;
    }
    return 0;
}
