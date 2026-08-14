#include <Trajectory/CartesianSpline.h>

#include <Eigen/Geometry>
#include <iostream>
#include <vector>

int main()
{
    using RobotLibrary::Model::Pose;
    using RobotLibrary::Trajectory::CartesianSpline;

    const Eigen::Quaterniond q0 = Eigen::Quaterniond::Identity();
    const Eigen::Quaterniond q1(Eigen::AngleAxisd(0.9, Eigen::Vector3d::UnitY()));
    const Eigen::Quaterniond q2(Eigen::AngleAxisd(1.4, Eigen::Vector3d(0.2, 0.9, -0.3).normalized()));
    Eigen::Quaterniond q1Negated = q1;
    q1Negated.coeffs() *= -1.0;

    const std::vector<double> times{0.0, 1.0, 2.0};
    const std::vector<Pose> canonical{
        {Eigen::Vector3d(0.0, 0.0, 0.0), q0},
        {Eigen::Vector3d(0.1, 0.2, 0.0), q1},
        {Eigen::Vector3d(0.2, 0.1, 0.1), q2},
    };
    const std::vector<Pose> signFlipped{
        canonical[0],
        {canonical[1].translation(), q1Negated},
        canonical[2],
    };

    CartesianSpline a(canonical, times, Eigen::Vector<double,6>::Zero());
    CartesianSpline b(signFlipped, times, Eigen::Vector<double,6>::Zero());
    for(const double time : {0.0, 0.25, 0.75, 1.0, 1.25, 1.75, 2.0})
    {
        const auto stateA = a.query_state(time);
        const auto stateB = b.query_state(time);
        const double rotationDifference =
            (stateA.pose.quaternion().toRotationMatrix()
             - stateB.pose.quaternion().toRotationMatrix()).norm();
        if(rotationDifference > 1e-10
           or (stateA.twist - stateB.twist).norm() > 1e-9
           or (stateA.acceleration - stateB.acceleration).norm() > 1e-8)
        {
            std::cerr << "CartesianSpline depends on quaternion sign at t=" << time
                      << ": rotation=" << rotationDifference
                      << ", twist=" << (stateA.twist - stateB.twist).norm()
                      << ", acceleration="
                      << (stateA.acceleration - stateB.acceleration).norm() << '\n';
            return 1;
        }
    }
    return 0;
}
