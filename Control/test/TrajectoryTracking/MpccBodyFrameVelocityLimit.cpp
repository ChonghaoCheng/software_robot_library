#include <Control/TrajectoryTracking/SerialLinkMPCC.h>

#include <Eigen/Geometry>
#include <array>
#include <cmath>
#include <iostream>

int main()
{
    constexpr double bound = 0.5;
    constexpr double tolerance = 1e-12;
    const Eigen::Matrix3d predictionRotation = Eigen::Matrix3d::Identity();
    const std::array<double,4> angles = {
        0.0, M_PI / 6.0, M_PI / 3.0, M_PI / 2.0};

    for(const double angle : angles)
    {
        const Eigen::Matrix3d stageRotation =
            Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitZ()).toRotationMatrix();
        const Eigen::Matrix<double,6,6> predictionToBody =
            RobotLibrary::Control::mpcc_prediction_twist_to_body_map(
                stageRotation, predictionRotation);

        Eigen::Vector<double,6> acceptedBody;
        acceptedBody << bound, -0.4, 0.3, -bound, 0.2, -0.1;
        const Eigen::Vector<double,6> acceptedPrediction =
            predictionToBody.transpose() * acceptedBody;
        const Eigen::Vector<double,6> recoveredAccepted =
            predictionToBody * acceptedPrediction;
        if((recoveredAccepted - acceptedBody).norm() > tolerance
           || (recoveredAccepted.array().abs() > bound + tolerance).any())
        {
            std::cerr << "accepted body twist changed at angle " << angle << '\n';
            return 1;
        }

        Eigen::Vector<double,6> rejectedBody = acceptedBody;
        rejectedBody(1) = bound + 1e-3;
        const Eigen::Vector<double,6> rejectedPrediction =
            predictionToBody.transpose() * rejectedBody;
        const Eigen::Vector<double,6> recoveredRejected =
            predictionToBody * rejectedPrediction;
        const bool directBodyAccepts =
            (rejectedBody.array().abs() <= bound + tolerance).all();
        const bool mappedMpccAccepts =
            (recoveredRejected.array().abs() <= bound + tolerance).all();
        if(directBodyAccepts != mappedMpccAccepts || mappedMpccAccepts)
        {
            std::cerr << "body-frame rejection mismatch at angle " << angle << '\n';
            return 1;
        }
    }

    std::cout << "body_frame_component_limit_m_per_s=" << bound << '\n'
              << "body_frame_component_limit_rad_per_s=" << bound << '\n'
              << "tested_rotations_deg=0,30,60,90\n";
    return 0;
}
