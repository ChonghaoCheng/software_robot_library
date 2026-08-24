/**
 * @file ContactGeometry.cpp
 * @brief Deterministic rigid-motion, SO(3)-branch, and impedance-error tests.
 */

#include <Control/Contact/SerialLinkMovingFrameMPC.h>
#include <ContactGeometry.h>

#include <Eigen/Geometry>

#include <cmath>
#include <iostream>
#include <string>

namespace {

int failures = 0;

void check(const bool condition, const std::string &message)
{
    if(!condition)
    {
        std::cerr << "FAIL: " << message << "\n";
        ++failures;
    }
}

bool near(const double lhs, const double rhs, const double tolerance = 1e-12)
{
    return std::abs(lhs - rhs) <= tolerance;
}

} // namespace

int main()
{
    using namespace RobotLibrary::Control;
    using namespace RobotLibrary::Control::detail;

    const ContactParameters defaults;
    const Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
    const Eigen::Vector3d defaultNormal =
        inward_contact_normal(identity, defaults.normalAxisInBoard);
    check((defaultNormal + Eigen::Vector3d::UnitY()).norm() <= 1e-15,
          "identity board maps the default inward normal to exactly -Y");

    const Eigen::Vector3d wristLinear(0.12, -0.03, 0.07);
    const Eigen::Vector3d wristAngular(0.0, 0.0, 2.0);
    const Eigen::Vector3d wristToPoint(0.0, 0.10, 0.0);
    const Eigen::Vector3d pointVelocity = offset_point_velocity(
        wristLinear, wristAngular, identity, wristToPoint);
    check((pointVelocity - Eigen::Vector3d(-0.08, -0.03, 0.07)).norm() <= 1e-15,
          "+Z angular velocity at a +Y offset contributes -X point velocity");

    Eigen::Vector<double,6> wristTwist;
    wristTwist << wristLinear, wristAngular;
    check((offset_point_velocity_map(identity, wristToPoint) * wristTwist
           - pointVelocity).norm() <= 1e-15,
          "[I, -skew(r)] maps the wrist-origin twist to the same point velocity");

    const Eigen::Matrix3d endpointRotation = Eigen::AngleAxisd(
        0.5 * M_PI, Eigen::Vector3d::UnitZ()).toRotationMatrix();
    const Eigen::Vector3d rotatedPointVelocity = offset_point_velocity(
        Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitY(),
        endpointRotation, wristToPoint);
    check((rotatedPointVelocity - Eigen::Vector3d(0.0, 0.0, 0.10)).norm() <= 2e-17,
          "the endpoint-fixed offset is rotated to world before applying omega cross r");

    const Eigen::Vector3d board0(0.4, -0.2, 0.7);
    const Eigen::Vector3d tool0 = board0 + Eigen::Vector3d(0.1, 0.02, 0.0);
    const double d0 = signed_normal_coordinate(defaultNormal, tool0, board0);
    const double movedForce = normal_force_from_coordinates(
        5.0, 1500.0, d0,
        signed_normal_coordinate(defaultNormal, tool0 + 0.001 * defaultNormal, board0));
    check(near(movedForce, 6.5, 1e-13),
          "positive tool displacement along default -Y increases compressive force");

    const Eigen::Vector3d commonTranslation(0.17, -0.31, 0.08);
    const double translatedCoordinate = signed_normal_coordinate(
        defaultNormal, tool0 + commonTranslation, board0 + commonTranslation);
    check(near(translatedCoordinate, d0, 2e-16),
          "common rigid translation leaves the signed normal coordinate invariant");
    check(near(normal_force_from_coordinates(5.0, 1500.0, d0, translatedCoordinate),
               5.0, 3e-13),
          "common rigid translation produces no artificial force variation");

    const Eigen::Vector3d rBoard(0.10, 0.02, 0.0);
    const Eigen::Matrix3d rotation = Eigen::AngleAxisd(
        5.0 * M_PI / 180.0, Eigen::Vector3d::UnitZ()).toRotationMatrix();
    const Eigen::Vector3d board1(-0.15, 0.22, 0.4);
    const Eigen::Vector3d toolAt0 = board0 + rBoard;
    const Eigen::Vector3d toolAt1 = board1 + rotation * rBoard;
    const Eigen::Vector3d normal0 = inward_contact_normal(identity, defaults.normalAxisInBoard);
    const Eigen::Vector3d normal1 = inward_contact_normal(rotation, defaults.normalAxisInBoard);
    const double rigidD0 = signed_normal_coordinate(normal0, toolAt0, board0);
    const double rigidD1 = signed_normal_coordinate(normal1, toolAt1, board1);
    check(near(rigidD1, rigidD0, 3e-17),
          "fixed board-frame contact offset is invariant under a five-degree board rotation");

    Eigen::Matrix3d positionMap = Eigen::Matrix3d::Identity();
    const AffineNormalForce affine = affine_normal_force(
        5.0, 1500.0, rigidD0, normal1, toolAt0, board1, positionMap);
    const double affineForce = affine.constant + affine.map.dot(toolAt1 - toolAt0);
    check(near(affineForce, 5.0, 5e-13),
          "the exact affine QP force expression is rigid-rotation invariant");
    check((normal1 - normal0).norm() > 0.08,
          "board rotation changes the future contact normal stage-by-stage");

    const Eigen::Vector3d axis = Eigen::Vector3d::UnitZ();
    const Eigen::Quaterniond current(Eigen::AngleAxisd(179.0 * M_PI / 180.0, axis));
    const Eigen::Quaterniond desired(Eigen::AngleAxisd(-179.0 * M_PI / 180.0, axis));
    const Eigen::Vector3d currentVector =
        RobotLibrary::Math::quaternion_to_rotation_vector(current);
    const Eigen::Vector3d localReference =
        local_orientation_reference(current, currentVector, desired);
    check(near((localReference - currentVector).norm(), 2.0 * M_PI / 180.0, 1e-12),
          "local orientation reference crosses the pi branch by two degrees, not 358");

    const RobotLibrary::Model::Pose origin(
        Eigen::Vector3d::Zero(), Eigen::Quaterniond::Identity());
    for(const double degrees : {10.0, 90.0})
    {
        const double radians = degrees * M_PI / 180.0;
        const RobotLibrary::Model::Pose rotated(
            Eigen::Vector3d::Zero(),
            Eigen::Quaterniond(Eigen::AngleAxisd(radians, Eigen::Vector3d::UnitX())));
        const Eigen::Vector<double,6> error = impedance_pose_error(origin, rotated);
        check(near(error.tail<3>().norm(), radians, 1e-12),
              "impedance orientation displacement has radian-valued angle-axis magnitude at "
              + std::to_string(static_cast<int>(degrees)) + " degrees");
        check(std::abs(error.tail<3>().norm() - std::sin(0.5 * radians)) > 1e-3,
              "impedance orientation displacement is not sin(theta/2)");
    }

    const RobotLibrary::Model::Pose movingSample(
        Eigen::Vector3d::Zero(), Eigen::Quaterniond::Identity());
    const RobotLibrary::Model::Pose propagated = propagate_rigid_frame(
        movingSample, Eigen::Vector3d(0.1, 0.0, 0.0), Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Zero(), 0.05);
    check(near(propagated.translation().x(), 0.005, 1e-15),
          "a 50 ms old 0.1 m/s sample is translated 0.005 m to control time");

    if(failures == 0)
    {
        std::cout << "contact_geometry_test PASS\n";
        return 0;
    }
    std::cerr << "contact_geometry_test FAIL (" << failures << " checks)\n";
    return 1;
}
