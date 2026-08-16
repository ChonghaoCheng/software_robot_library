#include <Control/TrajectoryTracking/SerialLinkCartesianMPCC.h>
#include <Math/MathFunctions.h>

#include <Eigen/Geometry>
#include <cmath>
#include <iostream>

namespace {

class CartesianMPCCProbe : public RobotLibrary::Control::SerialLinkCartesianMPCC
{
    public:
        using RobotLibrary::Control::SerialLinkCartesianMPCC::cartesian_path_tangent;
};

bool near(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const double tolerance)
{
    return (a - b).norm() <= tolerance;
}

} // namespace

int main()
{
    using RobotLibrary::Model::Pose;
    using RobotLibrary::Trajectory::CartesianSpline;

    const Pose start(Eigen::Vector3d::Zero(), Eigen::Quaterniond::Identity());
    const Pose finish(
        Eigen::Vector3d(0.4, -0.1, 0.2),
        Eigen::Quaterniond(Eigen::AngleAxisd(0.6, Eigen::Vector3d::UnitZ())));
    CartesianSpline trajectory(
        {start, finish}, {0.0, 2.0}, Eigen::Vector<double,6>::Zero());

    constexpr double progress = 0.4;
    constexpr double step = 1e-3;
    const Eigen::Vector<double,6> tangent =
        CartesianMPCCProbe::cartesian_path_tangent(trajectory, progress);
    const Eigen::Vector<double,6> intrinsic =
        trajectory.tangent_at_progress(progress, step);

    if(tangent.head<3>().norm() <= 1e-6)
    {
        std::cerr << "Cartesian position tangent is zero.\n";
        return 1;
    }
    if(tangent.tail<3>().norm() <= 1e-6)
    {
        std::cerr << "Rotating reference lost its angular feedforward.\n";
        return 2;
    }
    if(!near(tangent.tail<3>(), intrinsic.tail<3>(), 1e-10))
    {
        std::cerr << "Angular feedforward differs from the spline body tangent.\n";
        return 3;
    }

    return 0;
}
