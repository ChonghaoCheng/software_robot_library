#include <Model/Pose.h>

#include <Eigen/Geometry>

#include <cmath>
#include <iostream>
#include <random>
#include <string>

namespace
{
constexpr double kTolerance = 1e-11;

bool near(const Eigen::Vector3d &actual,
          const Eigen::Vector3d &expected,
          const double tolerance = kTolerance)
{
    return (actual - expected).norm() <= tolerance;
}

bool same_rotation(const Eigen::Quaterniond &actual,
                   const Eigen::Quaterniond &expected,
                   const double tolerance = kTolerance)
{
    return (actual.toRotationMatrix() - expected.toRotationMatrix()).norm() <= tolerance;
}

bool check_pose(const RobotLibrary::Model::Pose &actual,
                const RobotLibrary::Model::Pose &expected,
                const std::string &label)
{
    if(!near(actual.translation(), expected.translation()) ||
       !same_rotation(actual.quaternion(), expected.quaternion()))
    {
        std::cerr << label << " failed: translation error = "
                  << (actual.translation() - expected.translation()).norm()
                  << ", rotation-matrix error = "
                  << (actual.quaternion().toRotationMatrix() -
                      expected.quaternion().toRotationMatrix()).norm()
                  << '\n';
        return false;
    }
    return true;
}

bool check_inverse(RobotLibrary::Model::Pose transform,
                   const Eigen::Vector3d &point,
                   const std::string &label)
{
    RobotLibrary::Model::Pose inverse = transform.inverse();
    const RobotLibrary::Model::Pose identity;

    if(!check_pose(inverse * transform, identity, label + " inverse*T") ||
       !check_pose(transform * inverse, identity, label + " T*inverse"))
        return false;

    RobotLibrary::Model::Pose double_inverse = inverse.inverse();
    if(!check_pose(double_inverse, transform, label + " double inverse"))
        return false;

    const Eigen::Vector3d transformed_point = transform * point;
    const Eigen::Vector3d recovered_point = inverse * transformed_point;
    if(!near(recovered_point, point))
    {
        std::cerr << label << " point round trip failed: error = "
                  << (recovered_point - point).norm() << '\n';
        return false;
    }
    return true;
}
}

int main()
{
    const Eigen::Vector3d non_unit_axis(2.0, -3.0, 4.0);
    RobotLibrary::Model::Pose fixed_transform(
        Eigen::Vector3d(1.25, -2.5, 0.75),
        Eigen::Quaterniond(Eigen::AngleAxisd(1.7, non_unit_axis.normalized())));
    if(!check_inverse(fixed_transform, Eigen::Vector3d(-0.3, 2.1, 4.2), "fixed"))
        return 1;

    std::mt19937_64 generator(0x5E3A124ULL);
    std::uniform_real_distribution<double> component(-4.0, 4.0);
    std::uniform_real_distribution<double> angle_magnitude(0.2, 3.0);
    std::bernoulli_distribution angle_sign;

    for(int sample = 0; sample < 128; ++sample)
    {
        Eigen::Vector3d axis(component(generator), component(generator), component(generator));
        if(axis.squaredNorm() < 1e-6)
            axis = Eigen::Vector3d(1.0, 2.0, -1.0);

        const double angle = (angle_sign(generator) ? 1.0 : -1.0) *
                             angle_magnitude(generator);
        RobotLibrary::Model::Pose transform(
            Eigen::Vector3d(component(generator), component(generator), component(generator)),
            Eigen::Quaterniond(Eigen::AngleAxisd(angle, axis.normalized())));
        const Eigen::Vector3d point(component(generator), component(generator), component(generator));

        if(!check_inverse(transform, point, "random sample " + std::to_string(sample)))
            return 1;
    }

    return 0;
}
