#include <Control/TrajectoryTracking/SerialLinkMPCC.h>
#include <Trajectory/CartesianSpline.h>

#include <Eigen/Geometry>
#include <cmath>
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


    // FINAL-HARNESS-FREEZE-01: when the parent rotates, each horizon stage
    // must be expressed using the predicted parent orientation at that stage.
    const Eigen::Matrix4d parent0 =
        (Eigen::Translation3d(0.1, -0.2, 0.3)
         * Eigen::AngleAxisd(0.17, Eigen::Vector3d(1.0, 2.0, -1.0).normalized()))
            .matrix();
    Eigen::Vector<double,6> parentTwist;
    parentTwist << 0.01, -0.02, 0.03, 0.31, -0.23, 0.19;
    constexpr double dt = 0.002;
    double maxStageRotationError = 0.0;
    double maxMappedMotionError = 0.0;
    const Eigen::Matrix3d commonPredictionRotation =
        parent0.block<3,3>(0,0) * predictionRotation;
    for(int stage = 0; stage < 20; ++stage)
    {
        const double s = predictionFrameProgress
            + (stage + 0.5) * (stageProgress - predictionFrameProgress) / 20.0;
        const Eigen::Matrix4d pathPose = trajectory.pose_at_progress(s).as_matrix();
        const Eigen::Matrix4d predictedParent = parent0
            * RobotLibrary::Math::se3_exponential(stage * dt * parentTwist);
        const Eigen::Matrix3d expectedStageRotation =
            predictedParent.block<3,3>(0,0) * pathPose.block<3,3>(0,0);
        const Eigen::Matrix3d implementedStageRotation =
            RobotLibrary::Control::mpcc_stage_reference_rotation(
                predictedParent, pathPose);
        const double rotationError = Eigen::AngleAxisd(
            expectedStageRotation.transpose() * implementedStageRotation).angle();
        maxStageRotationError = std::max(maxStageRotationError,
                                         std::abs(rotationError));

        const Eigen::Vector<double,6> motion =
            trajectory.tangent_at_progress(s);
        Eigen::Vector<double,6> expectedMapped;
        const Eigen::Matrix3d stageToPrediction =
            commonPredictionRotation.transpose() * expectedStageRotation;
        expectedMapped.head<3>() = stageToPrediction * motion.head<3>();
        expectedMapped.tail<3>() = stageToPrediction * motion.tail<3>();
        const Eigen::Vector<double,6> implementedMapped =
            RobotLibrary::Control::mpcc_express_body_tangent_in_prediction_frame(
                motion, implementedStageRotation, commonPredictionRotation);
        maxMappedMotionError = std::max(
            maxMappedMotionError, (expectedMapped - implementedMapped).norm());
    }
    if(maxStageRotationError >= 1e-12 || maxMappedMotionError >= 1e-10)
    {
        std::cerr << "MPCC predicted-parent stage rotation mismatch: rotation="
                  << maxStageRotationError << " mapped_motion="
                  << maxMappedMotionError << '\n';
        return 1;
    }
    std::cout << "max_stage_rotation_log_error=" << maxStageRotationError << '\n'
              << "max_mapped_motion_error=" << maxMappedMotionError << '\n';
    return 0;
}
