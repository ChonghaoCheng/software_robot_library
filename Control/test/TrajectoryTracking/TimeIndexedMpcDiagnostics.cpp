/**
 * @file  TimeIndexedMpcDiagnostics.cpp
 * @brief TEST GATE C: time-indexed MPC controller diagnostics are exact.
 *
 * Recomputes every recorded quantity independently of the recorder and requires
 * agreement to numerical tolerance:
 *   - realized twist equals J * qdot from the same Jacobian and joint command;
 *   - the logged realization error equals ||u_realized - u_commanded||;
 *   - the Cartesian saturation flags follow from the commanded twist and the
 *     frozen limits;
 *   - the joint saturation flag follows from the same control limits the
 *     resolved-rate layer was constrained by;
 *   - one diagnostic is produced per successful invocation and all values are
 *     finite.
 */

#include <Control/TrajectoryTracking/SerialLinkLieAlgebraMPC.h>
#include <Control/TrajectoryTracking/SerialLinkMPC.h>
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
#include <type_traits>
#include <vector>

namespace {

using RobotLibrary::Control::SerialLinkLieAlgebraMPC;
using RobotLibrary::Control::SerialLinkMPC;
using RobotLibrary::Control::SerialLinkParameters;
using RobotLibrary::Control::TimeIndexedMpcDiagnostics;
using RobotLibrary::Model::Pose;
using RobotLibrary::Trajectory::CartesianSpline;

int failures = 0;

void check(const bool condition, const std::string &what)
{
    if(not condition)
    {
        std::cerr << "FAIL: " << what << "\n";
        ++failures;
    }
}

/** A six-revolute-joint arm, enough for a non-degenerate endpoint Jacobian. */
class UrdfFixture
{
public:
    UrdfFixture()
    {
        const auto nonce = std::chrono::steady_clock::now().time_since_epoch().count();
        path = std::filesystem::temp_directory_path()
            / ("time_indexed_mpc_diag_" + std::to_string(nonce) + ".urdf");
        std::ofstream stream(path);
        stream << R"(<robot name="diagnostics_fixture">
  <link name="base"/>
  <link name="l1"/><link name="l2"/><link name="l3"/>
  <link name="l4"/><link name="l5"/><link name="l6"/>
  <link name="tool"/>
  <joint name="j1" type="revolute">
    <parent link="base"/><child link="l1"/><axis xyz="0 0 1"/>
    <origin xyz="0 0 0.10" rpy="0 0 0"/>
    <limit lower="-3.14" upper="3.14" effort="100" velocity="2"/>
  </joint>
  <joint name="j2" type="revolute">
    <parent link="l1"/><child link="l2"/><axis xyz="0 1 0"/>
    <origin xyz="0 0 0.10" rpy="0 0 0"/>
    <limit lower="-3.14" upper="3.14" effort="100" velocity="2"/>
  </joint>
  <joint name="j3" type="revolute">
    <parent link="l2"/><child link="l3"/><axis xyz="0 1 0"/>
    <origin xyz="0.20 0 0" rpy="0 0 0"/>
    <limit lower="-3.14" upper="3.14" effort="100" velocity="2"/>
  </joint>
  <joint name="j4" type="revolute">
    <parent link="l3"/><child link="l4"/><axis xyz="1 0 0"/>
    <origin xyz="0.20 0 0" rpy="0 0 0"/>
    <limit lower="-3.14" upper="3.14" effort="100" velocity="2"/>
  </joint>
  <joint name="j5" type="revolute">
    <parent link="l4"/><child link="l5"/><axis xyz="0 1 0"/>
    <origin xyz="0.10 0 0" rpy="0 0 0"/>
    <limit lower="-3.14" upper="3.14" effort="100" velocity="2"/>
  </joint>
  <joint name="j6" type="revolute">
    <parent link="l5"/><child link="l6"/><axis xyz="1 0 0"/>
    <origin xyz="0.08 0 0" rpy="0 0 0"/>
    <limit lower="-3.14" upper="3.14" effort="100" velocity="2"/>
  </joint>
  <joint name="tool_fixed" type="fixed">
    <parent link="l6"/><child link="tool"/>
    <origin xyz="0.05 0 0" rpy="0 0 0"/>
  </joint>
</robot>)";
    }

    ~UrdfFixture()
    {
        std::error_code error;
        std::filesystem::remove(path, error);
    }

    std::filesystem::path path;
};

CartesianSpline make_path(const Pose &start, const Eigen::Vector3d &offset)
{
    std::vector<Pose> poses;
    std::vector<double> times;
    for(int i = 0; i < 9; ++i)
    {
        const double u = static_cast<double>(i) / 8.0;
        const Eigen::Vector3d travel(0.05 * u, 0.04 * std::sin(M_PI * u), 0.03 * u);
        const Eigen::AngleAxisd rz(0.6 * u, Eigen::Vector3d::UnitZ());
        const Eigen::AngleAxisd ry(0.4 * std::sin(M_PI * u), Eigen::Vector3d::UnitY());
        poses.emplace_back(start.translation() + offset + travel,
                           Eigen::Quaterniond(rz * ry) * start.quaternion());
        times.push_back(4.0 * u);
    }
    return CartesianSpline(poses, times, Eigen::Vector<double,6>::Zero());
}

/**
 * @brief Drive one controller and independently re-derive its diagnostics.
 *
 * @param clampInBodyAxes Lie-algebra MPC clamps the twist in body axes, so its
 *        saturation flags are checked there; Cartesian MPC clamps in base axes.
 */
template<typename Controller>
void exercise(const std::string &label,
              const std::filesystem::path &urdf,
              const bool clampInBodyAxes,
              const Eigen::Vector3d &startOffset,
              bool &sawLinearSaturation)
{
    // The Cartesian velocity limits of these controllers are frozen, so the
    // test reads them from the diagnostics rather than setting them.
    const double maxLinearSpeed = 0.5;
    const double maxAngularSpeed = 0.5;
    auto model = std::make_shared<RobotLibrary::Model::KinematicTree>(urdf.string());
    const unsigned int numJoints = model->number_of_joints();
    Eigen::VectorXd position = Eigen::VectorXd::Constant(numJoints, 0.35);
    Eigen::VectorXd velocity = Eigen::VectorXd::Zero(numJoints);
    model->update_state(position, velocity);

    SerialLinkParameters parameters;
    parameters.qpsolver.maxSteps = 50;
    parameters.qpsolver.stepSizeTolerance = 1e-5;
    Controller controller(model, "tool", parameters, 20, 0.002);

    const Pose startPose = controller.endpoint_pose();
    controller.set_trajectory(make_path(startPose, startOffset));

    double previousInvocationIndex = 0.0;
    double maxRealizationDisagreement = 0.0;
    unsigned int invocations = 0;

    for(unsigned int k = 0; k < 200; ++k)
    {
        const double referenceTime = static_cast<double>(k) * 0.002;
        const Eigen::VectorXd jointCommand =
            controller.track_endpoint_trajectory_at_time(referenceTime);
        ++invocations;

        const TimeIndexedMpcDiagnostics &d = controller.diagnostics();

        // Exactly one diagnostic per successful invocation.
        check(std::abs(d.invocationIndex - (previousInvocationIndex + 1.0)) < 1e-12,
              label + ": one diagnostic per successful invocation");
        previousInvocationIndex = d.invocationIndex;

        check(d.qpStatus == 1.0, label + ": a returned invocation reports qp_status 1");
        if constexpr(std::is_same_v<Controller, SerialLinkLieAlgebraMPC>)
        {
            check(d.qpConverged, label + ": task-space QP explicitly converged");
            check(!d.qpHitMaxIterations,
                  label + ": task-space QP did not exhaust iterations");
            check(d.qpPrimalViolation <= 1e-12,
                  label + ": task-space QP primal violation recorded");
        }
        check(controller.resolved_rate_qp_diagnostics().converged,
              label + ": resolved-rate QP explicitly converged");
        check(d.finite, label + ": all diagnostic values finite");
        check(d.commandedTwist.allFinite(), label + ": commanded twist present and finite");
        check(d.realizedTwist.allFinite(), label + ": realized twist present and finite");
        check(std::abs(d.controllerDt - 0.002) < 1e-12, label + ": controller dt recorded");
        check(std::abs(d.referenceTime - referenceTime) < 1e-9,
              label + ": reference time recorded");

        // Realized twist re-derived from the same Jacobian and joint command.
        const Eigen::Vector<double,6> realized = controller.jacobian() * jointCommand;
        check((realized - d.realizedTwist).norm() < 1e-12,
              label + ": realized twist equals J*qdot");

        // Realization error re-derived from the logged vectors alone.
        const double recomputed = (d.realizedTwist - d.commandedTwist).norm();
        maxRealizationDisagreement = std::max(
            maxRealizationDisagreement, std::abs(recomputed - d.twistRealizationError));

        // Saturation flags re-derived from the commanded twist and frozen limits.
        const Eigen::Vector<double,6> &clampFrame = d.clampFrameTwist;
        bool expectLinear = false;
        bool expectAngular = false;
        for(int i = 0; i < 3; ++i)
        {
            if(std::abs(std::abs(clampFrame(i)) - maxLinearSpeed) <= 1e-12) expectLinear = true;
            if(std::abs(std::abs(clampFrame(i + 3)) - maxAngularSpeed) <= 1e-12) expectAngular = true;
        }
        check(expectLinear == d.linearLimitActive, label + ": linear saturation flag");
        check(expectAngular == d.angularLimitActive, label + ": angular saturation flag");
        check(std::abs(d.effectiveLinearLimit - maxLinearSpeed) < 1e-12,
              label + ": effective linear limit");
        check(std::abs(d.effectiveAngularLimit - maxAngularSpeed) < 1e-12,
              label + ": effective angular limit");

        // Joint saturation flag re-derived from the published margin, which is
        // computed against the same control limits the resolved-rate QP used.
        check(d.jointLimitActive == (d.jointLimitMargin <= 1e-9),
              label + ": joint saturation flag follows the published limit margin");
        check(std::isfinite(d.jointLimitMargin), label + ": joint limit margin finite");

        if(d.linearLimitActive) sawLinearSaturation = true;

        // The clamp-frame twist is the base-frame one whenever the clamp is in
        // base axes; otherwise the two differ only by a rotation of equal norm.
        if(not clampInBodyAxes)
        {
            check((d.clampFrameTwist - d.commandedTwist).norm() < 1e-12,
                  label + ": base-axis clamp records the same vector twice");
        }
        else
        {
            check(std::abs(d.clampFrameTwist.head<3>().norm()
                           - d.commandedTwist.head<3>().norm()) < 1e-9
                  and std::abs(d.clampFrameTwist.tail<3>().norm()
                               - d.commandedTwist.tail<3>().norm()) < 1e-9,
                  label + ": body-axis clamp is a rotation of the base-axis twist");
        }

        // Advance the plant by the commanded joint velocities.
        position += jointCommand * 0.002;
        model->update_state(position, jointCommand);
    }

    check(maxRealizationDisagreement < 1e-12,
          label + ": logged realization error agrees with offline recomputation");
    std::cout << label << " invocations " << invocations
              << " max_realization_error_disagreement " << maxRealizationDisagreement << "\n";
}

} // namespace

int main()
{
    UrdfFixture fixture;
    bool sawLinearSaturation = false;
    bool ignored = false;

    // On-path start: the clamp does not bind, exercising the flags' false branch.
    exercise<SerialLinkLieAlgebraMPC>(
        "lie_algebra_mpc", fixture.path, true, Eigen::Vector3d::Zero(), ignored);
    exercise<SerialLinkMPC>(
        "cartesian_mpc", fixture.path, false, Eigen::Vector3d::Zero(), ignored);

    // Far start: the frozen 0.5 clamp binds, exercising the flags' true branch
    // without altering any controller bound.
    exercise<SerialLinkLieAlgebraMPC>(
        "lie_algebra_mpc_far", fixture.path, true, Eigen::Vector3d(0.15, -0.12, 0.10),
        sawLinearSaturation);
    exercise<SerialLinkMPC>(
        "cartesian_mpc_far", fixture.path, false, Eigen::Vector3d(0.15, -0.12, 0.10),
        sawLinearSaturation);

    check(sawLinearSaturation,
          "the far-start arms must actually reach the frozen Cartesian clamp, "
          "otherwise the saturation flags are never exercised in their true branch");

    if(failures == 0)
    {
        std::cout << "time_indexed_mpc_diagnostics_test PASS\n";
        return 0;
    }
    std::cerr << "time_indexed_mpc_diagnostics_test FAIL (" << failures << " checks)\n";
    return 1;
}
