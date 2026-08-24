/**
 * @file MovingFrameContactMpc.cpp
 * @brief Integration tests for force direction, board timestamps, validation, and fallback diagnostics.
 */

#include <Control/Contact/SerialLinkMovingFrameMPC.h>
#include <Model/KinematicTree.h>
#include <Trajectory/CartesianSpline.h>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using RobotLibrary::Control::ContactMode;
using RobotLibrary::Control::ContactParameters;
using RobotLibrary::Control::MovingFrameState;
using RobotLibrary::Control::SerialLinkMovingFrameMPC;
using RobotLibrary::Control::SerialLinkParameters;
using RobotLibrary::Model::Pose;
using RobotLibrary::Trajectory::CartesianSpline;

int failures = 0;

void check(const bool condition, const std::string &message)
{
    if(!condition)
    {
        std::cerr << "FAIL: " << message << "\n";
        ++failures;
    }
}

class UrdfFixture
{
public:
    UrdfFixture()
    {
        const auto nonce = std::chrono::steady_clock::now().time_since_epoch().count();
        path = std::filesystem::temp_directory_path()
            / ("moving_frame_contact_" + std::to_string(nonce) + ".urdf");
        std::ofstream stream(path);
        stream << R"(<robot name="contact_fixture">
  <link name="base"/>
  <link name="l1"/><link name="l2"/><link name="l3"/>
  <link name="l4"/><link name="l5"/><link name="l6"/>
  <link name="tool"/>
  <joint name="j1" type="revolute"><parent link="base"/><child link="l1"/><axis xyz="0 0 1"/><origin xyz="0 0 0.10" rpy="0 0 0"/><limit lower="-3.14" upper="3.14" effort="100" velocity="2"/></joint>
  <joint name="j2" type="revolute"><parent link="l1"/><child link="l2"/><axis xyz="0 1 0"/><origin xyz="0 0 0.10" rpy="0 0 0"/><limit lower="-3.14" upper="3.14" effort="100" velocity="2"/></joint>
  <joint name="j3" type="revolute"><parent link="l2"/><child link="l3"/><axis xyz="0 1 0"/><origin xyz="0.20 0 0" rpy="0 0 0"/><limit lower="-3.14" upper="3.14" effort="100" velocity="2"/></joint>
  <joint name="j4" type="revolute"><parent link="l3"/><child link="l4"/><axis xyz="1 0 0"/><origin xyz="0.20 0 0" rpy="0 0 0"/><limit lower="-3.14" upper="3.14" effort="100" velocity="2"/></joint>
  <joint name="j5" type="revolute"><parent link="l4"/><child link="l5"/><axis xyz="0 1 0"/><origin xyz="0.10 0 0" rpy="0 0 0"/><limit lower="-3.14" upper="3.14" effort="100" velocity="2"/></joint>
  <joint name="j6" type="revolute"><parent link="l5"/><child link="l6"/><axis xyz="1 0 0"/><origin xyz="0.08 0 0" rpy="0 0 0"/><limit lower="-3.14" upper="3.14" effort="100" velocity="2"/></joint>
  <joint name="tool_fixed" type="fixed"><parent link="l6"/><child link="tool"/><origin xyz="0.05 0 0" rpy="0 0 0"/></joint>
</robot>)";
    }

    ~UrdfFixture()
    {
        std::error_code error;
        std::filesystem::remove(path, error);
    }

    std::filesystem::path path;
};

std::shared_ptr<RobotLibrary::Model::KinematicTree>
make_model(const std::filesystem::path &urdf)
{
    auto model = std::make_shared<RobotLibrary::Model::KinematicTree>(urdf.string());
    model->update_state(Eigen::VectorXd::Constant(model->number_of_joints(), 0.35),
                        Eigen::VectorXd::Zero(model->number_of_joints()));
    return model;
}

void set_stationary_relative_trajectory(SerialLinkMovingFrameMPC &controller)
{
    const Pose pose = controller.endpoint_pose();
    controller.set_trajectory(CartesianSpline(
        pose, pose, Eigen::Vector<double,6>::Zero(), 0.0, 1.0));
}

std::vector<MovingFrameState> stationary_board_prediction()
{
    MovingFrameState board;
    board.pose = Pose(Eigen::Vector3d::Zero(), Eigen::Quaterniond::Identity());
    return {board, board};
}

ContactParameters force_parameters(const ContactMode mode)
{
    ContactParameters parameters;
    parameters.mode = mode;
    parameters.forceResponseGain = 1500.0;
    parameters.targetForce = 5.0;
    parameters.forceTolerance = 0.0;
    parameters.forceWeight = 1.0;
    parameters.slackWeight = 1e6;
    parameters.maxForceSlack = 50.0;
    parameters.tangentPositionWeight = 0.0;
    parameters.normalPositionWeight = 0.0;
    parameters.orientationWeight = 0.0;
    parameters.velocityWeight = 1e-3;
    parameters.deltaVelocityWeight = 0.0;
    parameters.normalAxisInBoard = Eigen::Vector3d::UnitZ();
    parameters.tangentAxisInBoard = Eigen::Vector3d::UnitX();
    return parameters;
}

void check_force_direction(const UrdfFixture &fixture,
                           const ContactMode mode,
                           const double measuredForce,
                           const double expectedSign,
                           const std::string &label)
{
    auto model = make_model(fixture.path);
    SerialLinkMovingFrameMPC controller(model, "tool", SerialLinkParameters(), 1, 0.01);
    set_stationary_relative_trajectory(controller);
    controller.set_contact_parameters(force_parameters(mode));
    controller.update_measured_normal_force(measuredForce);
    (void)controller.track_moving_frame_trajectory_at_time(0.0, stationary_board_prediction());

    const auto &diagnostics = controller.contact_diagnostics();
    check(diagnostics.solverSucceeded && !diagnostics.fallbackUsed,
          label + ": contact QP succeeds without fallback");
    check(expectedSign * diagnostics.firstCommand.head<3>().z() > 1e-6,
          label + ": first MPC velocity has the required inward-normal sign");
    check(expectedSign * (diagnostics.predictedFirstStepForce - measuredForce) > 0.0,
          label + ": first-step predicted force moves toward the target");
}

} // namespace

int main()
{
    const UrdfFixture fixture;

    check_force_direction(fixture, ContactMode::Loss, 3.0, 1.0, "Loss low force");
    check_force_direction(fixture, ContactMode::Loss, 7.0, -1.0, "Loss high force");
    check_force_direction(fixture, ContactMode::Constraint, 3.0, 1.0, "Constraint low force");
    check_force_direction(fixture, ContactMode::Constraint, 7.0, -1.0, "Constraint high force");

    {
        auto model = make_model(fixture.path);
        SerialLinkMovingFrameMPC controller(model, "tool", SerialLinkParameters(), 1, 0.01);
        set_stationary_relative_trajectory(controller);
        ContactParameters parameters = force_parameters(ContactMode::Loss);
        parameters.forceWeight = 0.0;
        controller.set_contact_parameters(parameters);
        controller.update_measured_normal_force(5.0);

        controller.update_board_pose(
            Pose(Eigen::Vector3d::Zero(), Eigen::Quaterniond::Identity()), 0.0);
        controller.update_board_pose(
            Pose(Eigen::Vector3d(0.005, 0.0, 0.0), Eigen::Quaterniond::Identity()), 0.05);
        controller.update_board_pose(
            Pose(Eigen::Vector3d(10.0, 0.0, 0.0), Eigen::Quaterniond::Identity()), 0.04);

        (void)controller.track_moving_frame_trajectory_at_time(0.10);
        const double expected = controller.endpoint_pose().translation().z();
        check(std::abs(controller.contact_diagnostics().currentSignedNormalCoordinate - expected) < 1e-12,
              "non-monotonic board sample is ignored without corrupting the predictor");

        // Repeat along +X so the 0.005 m age correction is directly observable.
        parameters.normalAxisInBoard = Eigen::Vector3d::UnitX();
        parameters.tangentAxisInBoard = Eigen::Vector3d::UnitY();
        controller.set_contact_parameters(parameters);
        (void)controller.track_moving_frame_trajectory_at_time(0.10);
        const double expectedNowCoordinate = controller.endpoint_pose().translation().x() - 0.010;
        check(std::abs(controller.contact_diagnostics().currentSignedNormalCoordinate
                       - expectedNowCoordinate) < 1e-12,
              "stage zero is propagated from the 50 ms old sample to the control timestamp");
    }

    {
        auto model = make_model(fixture.path);
        SerialLinkMovingFrameMPC controller(model, "tool", SerialLinkParameters(), 1, 0.01);
        ContactParameters invalid;
        invalid.forceResponseGain = 0.0;
        bool threw = false;
        try { controller.set_contact_parameters(invalid); }
        catch(const std::invalid_argument &) { threw = true; }
        check(threw, "zero contact stiffness is rejected");

        invalid = ContactParameters{};
        invalid.normalAxisInBoard.setZero();
        threw = false;
        try { controller.set_contact_parameters(invalid); }
        catch(const std::invalid_argument &) { threw = true; }
        check(threw, "zero board normal is rejected");
    }

    {
        auto model = make_model(fixture.path);
        SerialLinkMovingFrameMPC controller(model, "tool", SerialLinkParameters(), 1, 0.01);
        set_stationary_relative_trajectory(controller);
        ContactParameters impossible = force_parameters(ContactMode::Constraint);
        impossible.targetForce = 1e6;
        impossible.maxForceSlack = 0.0;
        controller.set_contact_parameters(impossible);
        controller.update_measured_normal_force(0.0);
        (void)controller.track_moving_frame_trajectory_at_time(0.0, stationary_board_prediction());
        check(!controller.contact_diagnostics().solverSucceeded
              && controller.contact_diagnostics().fallbackUsed,
              "infeasible contact QP exposes board-glued fallback in diagnostics");
    }

    if(failures == 0)
    {
        std::cout << "moving_frame_contact_mpc_test PASS\n";
        return 0;
    }
    std::cerr << "moving_frame_contact_mpc_test FAIL (" << failures << " checks)\n";
    return 1;
}
