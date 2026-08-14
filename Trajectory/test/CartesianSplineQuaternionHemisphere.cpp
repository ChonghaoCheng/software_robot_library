#include <Trajectory/CartesianSpline.h>

#include <Eigen/Geometry>
#include <iostream>
#include <string>
#include <vector>

namespace {

Eigen::Quaterniond negated(Eigen::Quaterniond quaternion)
{
    quaternion.coeffs() *= -1.0;
    return quaternion;
}

bool invariant(const std::vector<RobotLibrary::Model::Pose> &reference,
               const std::vector<RobotLibrary::Model::Pose> &candidate,
               const std::vector<double> &times,
               const std::string &label)
{
    using RobotLibrary::Trajectory::CartesianSpline;

    CartesianSpline a(reference, times, Eigen::Vector<double,6>::Zero());
    CartesianSpline b(candidate, times, Eigen::Vector<double,6>::Zero());
    for(const double time : {0.0, 0.25, 0.75, 1.0, 1.25, 1.75, 2.0})
    {
        const auto stateA = a.query_state(time);
        const auto stateB = b.query_state(time);
        const double rotationDifference =
            (stateA.pose.quaternion().toRotationMatrix()
             - stateB.pose.quaternion().toRotationMatrix()).norm();
        const double twistDifference = (stateA.twist - stateB.twist).norm();
        const double accelerationDifference =
            (stateA.acceleration - stateB.acceleration).norm();
        if(not std::isfinite(rotationDifference)
           or not std::isfinite(twistDifference)
           or not std::isfinite(accelerationDifference)
           or rotationDifference >= 1e-10
           or twistDifference >= 1e-9
           or accelerationDifference >= 1e-8)
        {
            std::cerr << label << " depends on quaternion sign at t=" << time
                      << ": rotation=" << rotationDifference
                      << ", twist=" << twistDifference
                      << ", acceleration=" << accelerationDifference << '\n';
            return false;
        }
    }
    return true;
}

} // namespace

int main()
{
    using RobotLibrary::Model::Pose;

    const Eigen::Quaterniond q0 = Eigen::Quaterniond::Identity();
    const Eigen::Quaterniond q1(Eigen::AngleAxisd(0.9, Eigen::Vector3d::UnitY()));
    const Eigen::Quaterniond q2(Eigen::AngleAxisd(1.4, Eigen::Vector3d(0.2, 0.9, -0.3).normalized()));

    const std::vector<double> times{0.0, 1.0, 2.0};
    const std::vector<Pose> a{
        {Eigen::Vector3d(0.0, 0.0, 0.0), q0},
        {Eigen::Vector3d(0.1, 0.2, 0.0), q1},
        {Eigen::Vector3d(0.2, 0.1, 0.1), q2},
    };
    const std::vector<Pose> b{
        {a[0].translation(), negated(q0)},
        {a[1].translation(), negated(q1)},
        {a[2].translation(), negated(q2)},
    };
    const std::vector<Pose> c{
        {a[0].translation(), negated(q0)}, a[1], a[2]
    };

    return invariant(a, b, times, "globally flipped sequence")
           and invariant(a, c, times, "first-only flipped sequence") ? 0 : 1;
}
