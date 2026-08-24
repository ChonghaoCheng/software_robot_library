/** @file ContactAdaptiveMPCC.cpp
 * Deterministic CA-MPCC-B geometry, QP, bound, warm-start and motion tests.
 */

#include <Control/Contact/SerialLinkContactAdaptiveMPCC.h>
#include <Model/KinematicTree.h>
#include <Trajectory/CartesianSpline.h>

#include <Eigen/Geometry>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

namespace {

using namespace RobotLibrary::Control;
using RobotLibrary::Model::Pose;
using RobotLibrary::Trajectory::CartesianSpline;
using RobotLibrary::Trajectory::CartesianTrajectoryFrameState;

int failures = 0;

void check(const bool condition, const std::string &message)
{
    if(!condition)
    {
        std::cerr << "FAIL: " << message << '\n';
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
            / ("contact_adaptive_mpcc_" + std::to_string(nonce) + ".urdf");
        std::ofstream stream(path);
        stream << R"(<robot name="ca_mpcc_b_fixture">
  <link name="base"/><link name="l1"/><link name="l2"/><link name="l3"/>
  <link name="l4"/><link name="l5"/><link name="l6"/><link name="tool"/>
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
make_model(const std::filesystem::path &path)
{
    auto model = std::make_shared<RobotLibrary::Model::KinematicTree>(path.string());
    model->update_state(Eigen::VectorXd::Constant(model->number_of_joints(), 0.35),
                        Eigen::VectorXd::Zero(model->number_of_joints()));
    return model;
}

template<class Controller>
void set_stationary_path(Controller &controller)
{
    const Pose pose = controller.endpoint_pose();
    controller.set_trajectory(CartesianSpline(
        pose, pose, Eigen::Vector<double,6>::Zero(), 0.0, 1.0));
}

ContactAdaptiveMpccParameters active_parameters()
{
    ContactAdaptiveMpccParameters parameters;
    parameters.targetForce = 5.0;
    parameters.deltaGain = 1e-3;
    parameters.forceDeadband = 0.1;
    parameters.deltaMinimum = -5e-4;
    parameters.deltaMaximum = 5e-4;
    parameters.deltaRateMaximum = 2e-3;
    parameters.forceRateWeight = 1e6;
    parameters.deltaRateSmoothWeight = 1e-3;
    parameters.compressionDirectionInParent = Eigen::Vector3d::UnitZ();
    return parameters;
}

Eigen::Matrix4d transform(const Eigen::Matrix3d &rotation,
                          const Eigen::Vector3d &translation)
{
    Eigen::Matrix4d result = Eigen::Matrix4d::Identity();
    result.block<3,3>(0,0) = rotation;
    result.block<3,1>(0,3) = translation;
    return result;
}

} // namespace

int main()
{
    const UrdfFixture fixture;

    // A: tau_delta is [R_d^T n_B;0] in stage body coordinates and is
    // identical to R_ref^T R_B n_B. Board rotation appears exactly once.
    {
        const Eigen::Matrix3d board =
            Eigen::AngleAxisd(0.63, Eigen::Vector3d(1.0, 2.0, -1.0).normalized())
                .toRotationMatrix();
        const Eigen::Matrix3d path =
            Eigen::AngleAxisd(-0.41, Eigen::Vector3d(-1.0, 1.0, 2.0).normalized())
                .toRotationMatrix();
        const Eigen::Vector3d normal = Eigen::Vector3d(0.2, -0.9, 0.4).normalized();
        const Eigen::Matrix3d reference = board * path;
        const auto body = contact_delta_tangent_body(path, normal);
        const auto prediction = contact_delta_tangent_prediction_frame(
            reference, board, normal);
        check((body - prediction).norm() <= 1e-12,
              "tau_delta frame formula reduces to R_d^T n_B");
        check(body.tail<3>().norm() == 0.0,
              "tau_delta has no angular component");
    }

    // B-D: force sign and deadband convention.
    {
        auto parameters = active_parameters();
        check(contact_preferred_delta_rate(parameters, 5.0) == 0.0,
              "zero force error gives zero preferred delta rate");
        check(contact_preferred_delta_rate(parameters, 3.0) > 0.0,
              "low force prefers positive compression-coordinate rate");
        check(contact_preferred_delta_rate(parameters, 7.0) < 0.0,
              "high force prefers negative compression-coordinate rate");
        check(contact_force_deadband(0.05, 0.1) == 0.0,
              "deadband removes small force errors");
    }

    // L: the 8D formulation with disabled adaptation/cost numerically
    // recovers the original 7D Basic MPCC.
    {
        auto baseModel = make_model(fixture.path);
        auto disabledModel = make_model(fixture.path);
        SerialLinkMPCC base(baseModel, "tool", SerialLinkParameters{}, 3, 0.01);
        SerialLinkContactAdaptiveMPCC disabled(
            disabledModel, "tool", SerialLinkParameters{}, 3, 0.01);
        ContactAdaptiveMpccParameters parameters = active_parameters();
        parameters.adaptationEnabled = false;
        parameters.contactCostEnabled = false;
        disabled.set_contact_adaptive_parameters(parameters);
        set_stationary_path(base);
        set_stationary_path(disabled);
        const Eigen::VectorXd baseCommand = base.step(0.01);
        const Eigen::VectorXd disabledCommand = disabled.step(0.01);
        check((baseCommand - disabledCommand).norm() <= 1e-9,
              "disabled 8D solver recovers Basic MPCC");
        check(disabled.diagnostics().stageControlDimension == 8,
              "disabled CA-MPCC-B still uses an 8D stage layout");
    }

    // E-G, M: rate bounds, cumulative delta bounds, interleaved warm shift,
    // and large force weight.
    {
        auto model = make_model(fixture.path);
        SerialLinkContactAdaptiveMPCC controller(
            model, "tool", SerialLinkParameters{}, 4, 0.01);
        auto parameters = active_parameters();
        controller.set_contact_adaptive_parameters(parameters);
        controller.set_measured_normal_force(3.0);
        set_stationary_path(controller);
        controller.set_contact_coordinate(parameters.deltaMaximum - 5e-6);
        (void)controller.step(0.01);
        const auto &contact = controller.contact_adaptive_diagnostics();
        check(contact.optimizedDeltaRates.size() == 4,
              "one optimized delta rate exists per horizon stage");
        check(contact.optimizedDeltaRates.cwiseAbs().maxCoeff()
                  <= parameters.deltaRateMaximum + 1e-9,
              "all optimized delta rates satisfy the rate bound");
        check(contact.predictedDelta.maxCoeff() <= parameters.deltaMaximum + 1e-9
              && contact.predictedDelta.minCoeff() >= parameters.deltaMinimum - 1e-9,
              "all cumulative horizon delta states satisfy coordinate bounds");
        const auto &mpcc = controller.diagnostics();
        check(mpcc.optimalHorizon.size() == 32 && mpcc.shiftedWarmStart.size() == 32,
              "8D horizon and warm start contain all variables");
        bool shifted = true;
        for(int stage = 0; stage < 3; ++stage)
        {
            shifted = shifted
                && (mpcc.shiftedWarmStart.segment(stage * 8, 8)
                    - mpcc.optimalHorizon.segment((stage + 1) * 8, 8)).norm() <= 1e-12;
        }
        check(shifted, "warm start shifts xi, s_dot, and delta_dot together");

        set_stationary_path(controller); // resets delta to the feasible centre
        controller.set_measured_normal_force(3.0);
        (void)controller.step(0.01);
        const auto &largeWeight = controller.contact_adaptive_diagnostics();
        check(std::abs(largeWeight.optimizedDeltaRate
                       - largeWeight.preferredDeltaRate) <= 1e-5,
              "very large force weight tracks the feasible preferred rate");
    }

    // N: with a tiny force-rate weight the trajectory objective is permitted
    // to reject the preferred virtual rate.
    {
        auto model = make_model(fixture.path);
        SerialLinkContactAdaptiveMPCC controller(
            model, "tool", SerialLinkParameters{}, 3, 0.01);
        auto parameters = active_parameters();
        parameters.forceRateWeight = 1e-12;
        parameters.deltaRateSmoothWeight = 0.0;
        controller.set_contact_adaptive_parameters(parameters);
        controller.set_measured_normal_force(3.0);
        set_stationary_path(controller);
        (void)controller.step(0.01);
        const auto &diagnostics = controller.contact_adaptive_diagnostics();
        check(std::abs(diagnostics.optimizedDeltaRate)
                  < 0.5 * std::abs(diagnostics.preferredDeltaRate),
              "very small force weight permits sacrificing force-rate preference");
    }

    // H-K: stationary, translating and rotating parents all use the one causal
    // parent-motion predictor. Rotation changes the expressed tangent but does
    // not create a second compensation term.
    {
        auto model = make_model(fixture.path);
        SerialLinkContactAdaptiveMPCC controller(
            model, "tool", SerialLinkParameters{}, 3, 0.01);
        auto parameters = active_parameters();
        parameters.forceRateWeight = 1e5;
        controller.set_contact_adaptive_parameters(parameters);
        controller.set_measured_normal_force(5.0);
        set_stationary_path(controller);

        CartesianTrajectoryFrameState frame;
        frame.transformInBase = Eigen::Matrix4d::Identity();
        controller.set_trajectory_frame(frame, 0.0);
        (void)controller.step_at_time(0.0, 0.01);
        check(!controller.diagnostics().parentFrameMotionActive,
              "stationary board has no parent motion compensation");

        frame.transformInBase = transform(
            Eigen::Matrix3d::Identity(), Eigen::Vector3d(0.001, 0.0, 0.0));
        controller.set_trajectory_frame(frame, 0.1);
        (void)controller.step_at_time(0.1, 0.01);
        check(controller.diagnostics().parentFrameMotionActive,
              "translating board activates causal parent prediction");
        check(controller.diagnostics().parentFrameBodyTwist.head<3>().norm() > 0.0,
              "translating board produces parent linear twist");

        frame.transformInBase = transform(
            Eigen::AngleAxisd(0.01, Eigen::Vector3d::UnitY()).toRotationMatrix(),
            Eigen::Vector3d(0.001, 0.0, 0.0));
        controller.set_trajectory_frame(frame, 0.2);
        (void)controller.step_at_time(0.2, 0.01);
        check(controller.diagnostics().parentFrameBodyTwist.tail<3>().norm() > 0.0,
              "rotating board produces parent angular twist");
        check(std::abs(controller.contact_adaptive_diagnostics()
                           .firstDeltaTangent.head<3>().norm() - 1.0) <= 1e-12,
              "rotating-board tau_delta remains a unit coordinate tangent");

        const Eigen::Matrix3d boardRotation =
            frame.transformInBase.block<3,3>(0,0);
        const Eigen::Matrix3d pathRotation =
            controller.endpoint_pose().quaternion().toRotationMatrix();
        const auto once = contact_delta_tangent_prediction_frame(
            boardRotation * pathRotation, boardRotation,
            parameters.compressionDirectionInParent);
        const auto reduced = contact_delta_tangent_body(
            pathRotation, parameters.compressionDirectionInParent);
        check((once - reduced).norm() <= 1e-12,
              "board rotation is not applied twice in tau_delta");
    }

    if(failures == 0)
    {
        std::cout << "contact_adaptive_mpcc_test PASS\n";
        return 0;
    }
    std::cerr << "contact_adaptive_mpcc_test FAIL (" << failures << " checks)\n";
    return 1;
}
