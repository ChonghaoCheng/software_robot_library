#include <Trajectory/CartesianSpline.h>

#include <Eigen/Geometry>

#include <algorithm>
#include <cstdio>
#include <vector>

int main()
{
    using RobotLibrary::Model::Pose;
    using RobotLibrary::Trajectory::CartesianSpline;

    const std::vector<Pose> poses{
        Pose(Eigen::Vector3d(0.1, -0.2, 0.3),
             Eigen::Quaterniond(Eigen::AngleAxisd(0.2, Eigen::Vector3d::UnitX()))),
        Pose(Eigen::Vector3d(0.3, 0.1, 0.4),
             Eigen::Quaterniond(Eigen::AngleAxisd(0.8, Eigen::Vector3d(1,2,3).normalized()))),
        Pose(Eigen::Vector3d(-0.1, 0.4, 0.2),
             Eigen::Quaterniond(Eigen::AngleAxisd(1.4, Eigen::Vector3d(2,-1,1).normalized())))};
    CartesianSpline spline(
        poses, {0.0, 1.3, 3.0}, Eigen::Vector<double,6>::Zero());

    double maximumDifference = 0.0;
    for(int i = 0; i <= 1000; ++i)
    {
        const double time = -0.1 + 3.2 * static_cast<double>(i) / 1000.0;
        const Eigen::Matrix4d full = spline.query_state(time).pose.as_matrix();
        const Eigen::Matrix4d poseOnly = spline.query_pose(time).as_matrix();
        maximumDifference = std::max(
            maximumDifference, (full - poseOnly).cwiseAbs().maxCoeff());
    }
    if(maximumDifference > 1e-15)
    {
        std::fprintf(stderr, "pose-only maximum difference %.17g exceeds 1e-15\n",
                     maximumDifference);
        return 1;
    }
    return 0;
}
