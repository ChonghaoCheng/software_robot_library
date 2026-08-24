/** @file ContactMPCC.cpp
 * Deterministic compatibility, force-direction, constraint, and timestamp tests.
 */

#include <Control/Contact/SerialLinkAdmittanceMPCC.h>
#include <Control/Contact/SerialLinkContactMPCC.h>
#include <Model/KinematicTree.h>
#include <Trajectory/CartesianSpline.h>

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
            / ("contact_mpcc_" + std::to_string(nonce) + ".urdf");
        std::ofstream stream(path);
        stream << R"(<robot name="contact_mpcc_fixture">
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

ContactParameters contact_parameters(const ContactMode mode)
{
    ContactParameters parameters;
    parameters.mode = mode;
    parameters.normalAxisInBoard = Eigen::Vector3d::UnitZ();
    parameters.tangentAxisInBoard = Eigen::Vector3d::UnitX();
    parameters.targetForce = 5.0;
    parameters.forceResponseGain = 1500.0;
    parameters.forceTolerance = 0.0;
    parameters.forceWeight = 10.0;
    parameters.slackWeight = 1e5;
    parameters.maxForceSlack = 50.0;
    parameters.pathLagPositionWeight = 0.0;
    parameters.tangentPositionWeight = 0.0;
    parameters.normalPositionWeight = 0.0;
    return parameters;
}

void check_coupled_direction(const UrdfFixture &fixture,
                             const ContactMode mode,
                             const double measuredForce,
                             const double sign,
                             const std::string &label)
{
    auto model = make_model(fixture.path);
    SerialLinkContactMPCC controller(
        model, "tool", SerialLinkParameters{}, 1, 0.01);
    set_stationary_path(controller);
    controller.set_contact_parameters(contact_parameters(mode));
    ContactState state;
    state.inContact = true;
    state.normalBase = Eigen::Vector3d::UnitZ();
    state.measuredNormalForce = measuredForce;
    controller.set_contact_state(state);
    (void)controller.step(0.01);

    const ContactMpcDiagnostics &diagnostics = controller.contact_diagnostics();
    check(diagnostics.solverSucceeded, label + ": augmented QP succeeds");
    check(sign * diagnostics.firstCommand.z() > 1e-7,
          label + ": commanded normal velocity points toward target force");
    check(sign * (diagnostics.predictedFirstStepForce - measuredForce) > 0.0,
          label + ": predicted force moves toward target");
    check(std::isfinite(diagnostics.realizedFirstStepForce)
          && std::isfinite(diagnostics.twistRealizationError),
          label + ": realized-twist force diagnostics are finite");
}

} // namespace

int main()
{
    const UrdfFixture fixture;

    // Identity postprocessor is an exact compatibility gate for legacy MPCC.
    {
        auto baseModel = make_model(fixture.path);
        auto composedModel = make_model(fixture.path);
        SerialLinkMPCC base(baseModel, "tool", SerialLinkParameters{}, 2, 0.01);
        SerialLinkAdmittanceMPCC composed(
            composedModel, "tool", SerialLinkParameters{}, 2, 0.01);
        set_stationary_path(base);
        set_stationary_path(composed);
        ContactReference disabled;
        disabled.enabled = false;
        composed.set_contact_reference(disabled);
        const Eigen::VectorXd baseCommand = base.step(0.01);
        const Eigen::VectorXd composedCommand = composed.step(0.01);
        check((baseCommand - composedCommand).norm() <= 1e-12,
              "disabled admittance composition preserves legacy MPCC numerically");
    }

    {
        auto model = make_model(fixture.path);
        SerialLinkAdmittanceMPCC controller(
            model, "tool", SerialLinkParameters{}, 1, 0.01);
        set_stationary_path(controller);
        ContactState state;
        state.normalBase = Eigen::Vector3d::UnitZ();
        state.inContact = true;
        state.measuredNormalForce = 3.0;
        controller.set_contact_state(state);
        ContactReference reference;
        reference.targetNormalForce = 5.0;
        reference.forceResponseGain = 1000.0;
        reference.maxNormalVelocity = 1.0;
        controller.set_contact_reference(reference);
        (void)controller.step(0.01);
        check(std::abs(controller.diagnostics().commandedBaseTwist.z() - 0.002) <= 1e-12,
              "admittance comparator replaces MPCC normal velocity after the QP");
    }

    check_coupled_direction(fixture, ContactMode::Loss, 3.0, 1.0, "Loss low force");
    check_coupled_direction(fixture, ContactMode::Loss, 7.0, -1.0, "Loss high force");
    check_coupled_direction(fixture, ContactMode::Constraint, 3.0, 1.0, "Constraint low force");
    check_coupled_direction(fixture, ContactMode::Constraint, 7.0, -1.0, "Constraint high force");

    {
        auto model = make_model(fixture.path);
        SerialLinkContactMPCC controller(model, "tool", SerialLinkParameters{}, 1, 0.01);
        set_stationary_path(controller);
        ContactParameters parameters = contact_parameters(ContactMode::Loss);
        parameters.approachWhenNotInContact = true;
        controller.set_contact_parameters(parameters);
        ContactState state;
        state.inContact = false;
        state.normalBase = Eigen::Vector3d::UnitZ();
        state.measuredNormalForce = 0.0;
        controller.set_contact_state(state);
        (void)controller.step(0.01);
        check(controller.contact_diagnostics().contactModeActive,
              "explicit approach option keeps the force objective active before contact");
        check(controller.contact_diagnostics().firstCommand.z() > 1e-7,
              "pre-contact force objective commands motion toward the known surface");
    }

    {
        CausalParentFrameMotion motion;
        motion.update(Eigen::Matrix4d::Identity(), 0.0);
        Eigen::Matrix4d translated = Eigen::Matrix4d::Identity();
        translated(0,3) = 1.0;
        motion.update(translated, 1.0);
        const Eigen::Matrix4d atControlTime =
            motion.predicted_pose_with_age(0.5, 0, 0.1);
        check(std::abs(atControlTime(0,3) - 1.5) <= 1e-12,
              "parent predictor compensates sample age before horizon stage zero");
    }

    {
        auto model = make_model(fixture.path);
        SerialLinkContactMPCC controller(model, "tool", SerialLinkParameters{}, 1, 0.01);
        ContactParameters invalid;
        invalid.tangentAxisInBoard = invalid.normalAxisInBoard;
        bool threw = false;
        try { controller.set_contact_parameters(invalid); }
        catch(const std::invalid_argument &) { threw = true; }
        check(threw, "parallel contact normal/tangent axes are rejected");

        ContactState state;
        state.measuredNormalForce = -1.0;
        threw = false;
        try { controller.set_contact_state(state); }
        catch(const std::invalid_argument &) { threw = true; }
        check(threw, "negative compressive-force measurements are rejected consistently");
    }

    if(failures == 0)
    {
        std::cout << "contact_mpcc_test PASS\n";
        return 0;
    }
    std::cerr << "contact_mpcc_test FAIL (" << failures << " checks)\n";
    return 1;
}
