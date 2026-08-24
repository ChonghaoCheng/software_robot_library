/** @file PredictiveContactMPCC.cpp */

#include <Control/Contact/PredictiveContactModel.h>
#include <Control/Contact/SerialLinkPredictiveContactMPCC.h>
#include <Control/TrajectoryTracking/SerialLinkKinematic.h>
#include <Math/MathFunctions.h>
#include <Model/KinematicTree.h>
#include <Trajectory/CartesianSpline.h>

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <string>

namespace {

using namespace RobotLibrary::Control;
using RobotLibrary::Model::Pose;
using RobotLibrary::Trajectory::CartesianSpline;
using RobotLibrary::Trajectory::CartesianTrajectoryFrameState;

int failures = 0;
double maximumAffineForceError = 0.0;

void check(const bool condition, const std::string &message)
{
    if(!condition)
    {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

Eigen::Matrix4d transform(const Eigen::Matrix3d &rotation,
                          const Eigen::Vector3d &translation)
{
    Eigen::Matrix4d result = Eigen::Matrix4d::Identity();
    result.block<3,3>(0,0) = rotation;
    result.block<3,1>(0,3) = translation;
    return result;
}

PredictiveContactKinematics kinematics(
    const Eigen::Vector<double,6> &boardBodyTwist =
        Eigen::Vector<double,6>::Zero())
{
    PredictiveContactKinematics result;
    result.horizon = 12;
    result.dt = 0.002;
    result.expectedDt = 0.002;
    result.predictionRotation = Eigen::AngleAxisd(
        0.31, Eigen::Vector3d(1.0, -2.0, 0.7).normalized()).toRotationMatrix();
    result.endpointRotationBase = Eigen::AngleAxisd(
        -0.47, Eigen::Vector3d(-0.2, 1.0, 0.5).normalized()).toRotationMatrix();
    result.endpointPositionBase = Eigen::Vector3d(0.42, -0.18, 0.63);
    result.contactOffsetEndpoint =
        Eigen::Vector3d(0.0005265895, 0.0999391081, -0.0399956612);
    result.compressionDirectionParent = Eigen::Vector3d(0.0, -1.0, 0.0);
    const Eigen::Matrix4d initial = transform(
        Eigen::AngleAxisd(0.21, Eigen::Vector3d::UnitX()).toRotationMatrix(),
        Eigen::Vector3d(0.4, -0.25, 0.6));
    for(int stage = 0; stage <= result.horizon; ++stage)
    {
        result.parentTransforms.push_back(
            initial * RobotLibrary::Math::se3_exponential(
                static_cast<double>(stage) * result.dt * boardBodyTwist));
    }
    return result;
}

Eigen::VectorXd random_decision(const int size)
{
    std::mt19937 generator(1706);
    std::uniform_real_distribution<double> distribution(-0.04, 0.04);
    Eigen::VectorXd decision(size);
    for(int index = 0; index < size; ++index)
    {
        decision(index) = distribution(generator);
    }
    for(int stage = 0; stage < size / 7; ++stage)
    {
        decision(7 * stage + 6) = 0.5;
    }
    return decision;
}

void check_affine_case(const PredictiveContactKinematics &input,
                       const Eigen::VectorXd &decision,
                       const std::string &label)
{
    constexpr double force = 2.5;
    constexpr double stiffness = 21303.75539503847;
    const PredictiveContactAffineModel affine =
        build_predictive_contact_affine_model(input);
    const Eigen::VectorXd condensed = predict_contact_force_affine(
        affine, decision, force, stiffness);
    const PredictiveContactRollout explicitRollout =
        rollout_predictive_contact_explicit(
            input, decision, force, stiffness);
    maximumAffineForceError = std::max(
        maximumAffineForceError,
        (condensed - explicitRollout.force).cwiseAbs().maxCoeff());
    check((condensed - explicitRollout.force).cwiseAbs().maxCoeff() < 1e-10,
          label + ": condensed force equals explicit rollout near machine precision");
    check((affine.relativeVelocityMap * decision
           + affine.relativeVelocityOffset
           - explicitRollout.relativeNormalVelocity).cwiseAbs().maxCoeff() < 1e-13,
          label + ": condensed relative velocity equals explicit rollout");
    check((affine.commandedRobotNormalVelocityMap * decision
           - explicitRollout.commandedRobotNormalVelocity).cwiseAbs().maxCoeff() < 1e-13,
          label + ": commanded robot-normal velocity equals explicit rollout");
    check((affine.realizedRobotNormalVelocityMap * decision
           + affine.realizedRobotNormalVelocityOffset
           - explicitRollout.realizedRobotNormalVelocity).cwiseAbs().maxCoeff() < 1e-13,
          label + ": realized robot-normal velocity equals explicit rollout");
}

class UrdfFixture
{
public:
    UrdfFixture()
    {
        const auto nonce = std::chrono::steady_clock::now().time_since_epoch().count();
        path = std::filesystem::temp_directory_path()
            / ("predictive_contact_mpcc_" + std::to_string(nonce) + ".urdf");
        std::ofstream stream(path);
        stream << R"(<robot name="predictive_contact_fixture">
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
void set_path(Controller &controller)
{
    const Pose pose = controller.endpoint_pose();
    controller.set_trajectory(CartesianSpline(
        pose, pose, Eigen::Vector<double,6>::Zero(), 0.0, 1.0));
    CartesianTrajectoryFrameState frame;
    frame.transformInBase = Eigen::Matrix4d::Identity();
    controller.set_trajectory_frame(frame, 0.0);
}

template<class Controller>
void set_moving_path(Controller &controller)
{
    const Pose start = controller.endpoint_pose();
    const Pose end(start.translation() + Eigen::Vector3d(0.02, 0.0, 0.0),
                   start.quaternion());
    controller.set_trajectory(CartesianSpline(
        start, end, Eigen::Vector<double,6>::Zero(), 0.0, 2.0));
    CartesianTrajectoryFrameState frame;
    frame.transformInBase = Eigen::Matrix4d::Identity();
    controller.set_trajectory_frame(frame, 0.0);
}

} // namespace

int main()
{
    // 1-9: point map, frozen-orientation rollout, translation, pitch,
    // Omega*M_P coupling, no double subtraction, zero/random decisions and sign.
    const PredictiveContactKinematics stationary = kinematics();
    Eigen::VectorXd zero = Eigen::VectorXd::Zero(84);
    check_affine_case(stationary, zero, "stationary zero U");

    Eigen::VectorXd angular = Eigen::VectorXd::Zero(84);
    angular(3) = 0.03;
    angular(4) = -0.02;
    angular(5) = 0.04;
    check_affine_case(stationary, angular, "stationary angular wrist motion");

    Eigen::Vector<double,6> translationTwist = Eigen::Vector<double,6>::Zero();
    translationTwist.head<3>() = Eigen::Vector3d(0.0, -2.5e-5, 0.0);
    check_affine_case(kinematics(translationTwist), random_decision(84),
                      "board normal translation with random bounded U");

    Eigen::Vector<double,6> pitchTwist = Eigen::Vector<double,6>::Zero();
    pitchTwist.tail<3>() = Eigen::Vector3d(0.7e-3, 0.0, 0.0);
    const PredictiveContactKinematics pitch = kinematics(pitchTwist);
    Eigen::VectorXd tangential = random_decision(84);
    for(int stage = 0; stage < 12; ++stage)
    {
        tangential.segment<3>(7 * stage) = Eigen::Vector3d(0.02, 0.0, 0.0);
    }
    check_affine_case(pitch, tangential,
                      "simultaneous robot tangential motion and board pitch");

    const PredictiveContactAffineModel pitchMaps =
        build_predictive_contact_affine_model(pitch);
    const int stage = 7;
    Eigen::MatrixXd directVelocity = Eigen::MatrixXd::Zero(3, 84);
    directVelocity.block(0, 7 * stage, 3, 7) =
        pitchMaps.pointVelocityMaps[static_cast<size_t>(stage)];
    const Eigen::RowVectorXd omittedOmegaTerm =
        pitchMaps.normalsBase[static_cast<size_t>(stage)].transpose()
        * directVelocity;
    check((pitchMaps.relativeVelocityMap.row(stage) - omittedOmegaTerm).norm() > 1e-9,
          "board pitch includes the nonzero -n^T Omega_B M_P U term");

    const PredictiveContactAffineModel translatingMaps =
        build_predictive_contact_affine_model(kinematics(translationTwist));
    const double expectedStationaryRobotRelativeVelocity =
        -translatingMaps.normalsBase.front().dot(
            translatingMaps.boardLinearVelocitiesBase.front()
            + translatingMaps.boardAngularVelocitiesBase.front().cross(
                translatingMaps.initialContactPointBase
                - kinematics(translationTwist).parentTransforms.front().block<3,1>(0,3)));
    check(std::abs(translatingMaps.relativeVelocityOffset(0)
                   - expectedStationaryRobotRelativeVelocity) < 1e-14,
          "board surface velocity is subtracted exactly once");

    const PredictiveContactAffineModel stationaryMaps =
        build_predictive_contact_affine_model(stationary);
    Eigen::VectorXd compression = Eigen::VectorXd::Zero(84);
    const Eigen::Vector3d normal = stationaryMaps.normalsBase.front();
    compression.head<3>() = stationary.predictionRotation.transpose() * normal * 1e-3;
    const double positive = (stationaryMaps.relativeVelocityMap * compression)(0);
    compression.head<3>() *= -1.0;
    const double negative = (stationaryMaps.relativeVelocityMap * compression)(0);
    check(positive > 0.0 && negative < 0.0,
          "positive/negative relative normal commands preserve compression sign");

    const Eigen::Matrix<double,3,7> &pointMap = stationaryMaps.pointVelocityMaps.front();
    Eigen::Vector<double,7> input = Eigen::Vector<double,7>::Zero();
    input.head<3>() = Eigen::Vector3d(0.1, -0.03, 0.07);
    input.segment<3>(3) = Eigen::Vector3d(0.0, 0.0, 2.0);
    const Eigen::Vector3d originLinear =
        stationary.predictionRotation * input.head<3>();
    const Eigen::Vector3d angularBase =
        stationary.predictionRotation * input.segment<3>(3);
    const Eigen::Vector3d offsetBase =
        stationary.endpointRotationBase * stationary.contactOffsetEndpoint;
    check((pointMap * input
           - (originLinear + angularBase.cross(offsetBase))).norm() < 1e-14,
          "contact-point affine velocity map matches omega cross r helper convention");
    Eigen::Vector<double,6> measuredWristTwist;
    measuredWristTwist << originLinear, angularBase;
    check((wrist_twist_to_contact_point_velocity(
               measuredWristTwist, stationary.endpointRotationBase,
               stationary.contactOffsetEndpoint)
           - (originLinear + angularBase.cross(offsetBase))).norm() < 1e-14,
          "measured wrist twist reconstructs the fixed-offset contact-point velocity");

    // E07-A: exact scalar realization condensation for all candidate delay
    // branches, measured nonzero z0, known past commands, angular lever arm,
    // translation, pitch, zero U, and random bounded U.
    for(int delay = 0; delay <= 5; ++delay)
    {
        for(const auto &base : {stationary, kinematics(translationTwist), pitch})
        {
            PredictiveContactKinematics aware = base;
            aware.actuationAware = true;
            aware.realizationAutoregressive = 0.9809437967444158;
            aware.realizationInputGain = 0.007248865304869999;
            aware.realizationDelay = delay;
            aware.initialRealizedRobotNormalVelocity = 3.7e-5;
            for(int sample = 0; sample < delay; ++sample)
                aware.pastRobotNormalCommands.push_back(
                    -2.0e-5 + 7.0e-6 * static_cast<double>(sample));
            check_affine_case(aware, zero,
                              "actuation-aware zero U delay " + std::to_string(delay));
            check_affine_case(aware, random_decision(84),
                              "actuation-aware random U delay " + std::to_string(delay));
        }
    }

    // 10-12 and 17: normalized cost is symmetric/PSD with the exact gradient;
    // M2 and M3 coincide at the same stiffness.
    const Eigen::MatrixXd forceMap = 21303.75539503847
        * stationaryMaps.penetrationIncrementMap;
    const Eigen::VectorXd forceConstant = Eigen::VectorXd::Constant(12, 2.35)
        + 21303.75539503847 * stationaryMaps.penetrationIncrementOffset;
    Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(84, 84);
    Eigen::VectorXd gradient = Eigen::VectorXd::Zero(84);
    add_normalized_predictive_force_cost(
        forceMap, forceConstant, 2.5, 3.0, hessian, gradient);
    check((hessian - hessian.transpose()).norm() < 1e-14,
          "normalized force Hessian is symmetric");
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenvalues(hessian);
    check(eigenvalues.eigenvalues().minCoeff() >= -1e-10,
          "normalized force Hessian is PSD within numerical tolerance");
    const Eigen::VectorXd direction = random_decision(84);
    const double epsilon = 1e-7;
    const auto objective = [&](const Eigen::VectorXd &u)
    {
        const Eigen::VectorXd residual =
            (forceMap * u + forceConstant
             - Eigen::VectorXd::Constant(12, 2.5)) / 2.5;
        return 3.0 * residual.squaredNorm();
    };
    const double finiteDifference =
        (objective(epsilon * direction) - objective(-epsilon * direction))
        / (2.0 * epsilon);
    check(std::abs(finiteDifference - gradient.dot(direction)) < 1e-7,
          "normalized force-cost gradient matches directional finite difference");
    const Eigen::VectorXd velocity =
        stationaryMaps.relativeVelocityMap * random_decision(84)
        + stationaryMaps.relativeVelocityOffset;
    Eigen::VectorXd m2 = Eigen::VectorXd::Constant(12, 2.5);
    double force = 2.5;
    for(int index = 0; index < 12; ++index)
    {
        force += 21303.75539503847 * 0.002 * velocity(index);
        m2(index) = force;
    }
    check((m2 - predict_contact_force_affine(
               stationaryMaps, random_decision(84), 2.5,
               21303.75539503847)).norm() < 1e-12,
          "M2 and M3 are identical for the same stiffness and initialization");

    // 13-16 and 18: active soft bounds, disabled recovery, explicit stale
    // policy and loud dt mismatch.
    const UrdfFixture fixture;
    {
        auto model = make_model(fixture.path);
        SerialLinkParameters solverParameters;
        solverParameters.mpccQpStepSizeTolerance = 1e-4;
        solverParameters.qpsolver.maxSteps = 100;
        SerialLinkPredictiveContactMPCC controller(
            model, "tool", solverParameters, 3, 0.002);
        PredictiveContactMpccParameters parameters;
        parameters.mode = PredictiveContactMode::Active;
        parameters.compressionDirectionParent = Eigen::Vector3d::UnitZ();
        parameters.contactOffsetEndpoint.setZero();
        parameters.actuationAware = true;
        parameters.realizationDelay = 4;
        parameters.normalActionGuardEnabled = true;
        controller.set_predictive_contact_parameters(parameters);
        controller.set_normal_realization_state(
            0.0, std::vector<double>(4, 0.0));
        set_path(controller);
        controller.set_force_measurement(3.25, true);
        controller.set_normal_realization_state(
            0.0, std::vector<double>(4, 0.0));
        (void)controller.step_at_time(0.0, 0.002);
        const auto &diagnostics = controller.predictive_contact_diagnostics();
        check(diagnostics.predictedForce.size() == 3
              && diagnostics.predictedForce.maxCoeff()
                   <= parameters.maximumForce + diagnostics.maximumForceSlack + 1e-7
              && diagnostics.predictedForce.minCoeff()
                   >= parameters.minimumForce - diagnostics.maximumForceSlack - 1e-7,
              "soft force upper/lower affine constraints hold");
        check(diagnostics.predictedPenetrationIncrement.size() == 3
              && diagnostics.predictedPenetrationIncrement.cwiseAbs().maxCoeff()
                  <= parameters.maximumPenetrationIncrement
                     + diagnostics.maximumPenetrationSlack + 1e-9,
              "soft penetration trust-region constraints hold");
        check(controller.diagnostics().stageControlDimension == 7,
              "actuation-aware prediction does not change the stage decision dimension");
        check(diagnostics.predictedCommandedRobotNormalVelocity.size() == 3
              && diagnostics.predictedCommandedRobotNormalVelocity.cwiseAbs().maxCoeff()
                  <= parameters.maximumRobotNormalCommand + 1e-10,
              "robot-normal command magnitude constraints hold");
        Eigen::VectorXd slew = diagnostics.predictedCommandedRobotNormalVelocity;
        for(int index = slew.size() - 1; index > 0; --index)
            slew(index) -= diagnostics.predictedCommandedRobotNormalVelocity(index - 1);
        check(slew.size() == 3 && slew.cwiseAbs().maxCoeff()
                  <= parameters.maximumRobotNormalCommandStep + 1e-10,
              "robot-normal command slew constraints hold");
    }

    {
        auto baseModel = make_model(fixture.path);
        auto disabledModel = make_model(fixture.path);
        SerialLinkMPCC base(baseModel, "tool", SerialLinkParameters{}, 3, 0.002);
        SerialLinkPredictiveContactMPCC disabled(
            disabledModel, "tool", SerialLinkParameters{}, 3, 0.002);
        PredictiveContactMpccParameters parameters;
        parameters.mode = PredictiveContactMode::Disabled;
        disabled.set_predictive_contact_parameters(parameters);
        set_path(base);
        set_path(disabled);
        const Eigen::VectorXd baseCommand = base.step_at_time(0.0, 0.002);
        const Eigen::VectorXd disabledCommand = disabled.step_at_time(0.0, 0.002);
        check((baseCommand - disabledCommand).norm() <= 1e-12,
              "disabled predictive force cost recovers Basic MPCC numerically");
    }

    {
        auto model = make_model(fixture.path);
        SerialLinkPredictiveContactMPCC controller(
            model, "tool", SerialLinkParameters{}, 3, 0.002);
        PredictiveContactMpccParameters parameters;
        parameters.mode = PredictiveContactMode::Shadow;
        parameters.compressionDirectionParent = Eigen::Vector3d::UnitZ();
        parameters.contactOffsetEndpoint.setZero();
        controller.set_predictive_contact_parameters(parameters);
        set_path(controller);
        controller.set_force_measurement(2.5, true);
        (void)controller.step_at_time(0.0, 0.002);
        check(controller.predictive_contact_diagnostics().predictedForce.size() == 3,
              "valid force anchors the shadow prediction");
        controller.set_force_measurement(2.5, false);
        (void)controller.step_at_time(0.002, 0.002);
        check(not controller.predictive_contact_diagnostics().forceValid
              && controller.predictive_contact_diagnostics().predictedForce.size() == 0,
              "stale force explicitly disables prediction instead of reusing the old anchor");

        controller.set_force_measurement(0.1, true);
        (void)controller.step_at_time(0.003, 0.002);
        check(controller.predictive_contact_diagnostics().forceValid
              && not controller.predictive_contact_diagnostics().contactModelValid
              && controller.predictive_contact_diagnostics().predictedForce.size() == 0,
              "fresh near-zero force disables M3 instead of declaring valid contact");

        controller.set_force_measurement(2.5, true);
        bool threw = false;
        try { (void)controller.step_at_time(0.004, 0.003); }
        catch(const std::invalid_argument &) { threw = true; }
        check(threw, "identified/control dt mismatch fails loudly");
    }

    // E06-Q: the legacy absolute tolerance suppresses small resolved-rate
    // commands, while a layer-specific tolerance recovers the high-accuracy
    // solution without violating constraints or destabilizing the damped branch.
    {
        Eigen::Vector<double,6> request = Eigen::Vector<double,6>::Zero();
        request(0) = 5e-5;
        SerialLinkParameters legacyParameters;
        legacyParameters.qpsolver.stepSizeTolerance = 5e-2;
        legacyParameters.qpsolver.maxSteps = 100;
        auto legacyModel = make_model(fixture.path);
        RobotLibrary::Control::SerialLinkKinematic legacy(
            legacyModel, "tool", legacyParameters);
        const Eigen::VectorXd legacyCommand = legacy.resolve_endpoint_twist(request);
        check(legacyCommand.norm() == 0.0,
              "legacy 0.05 tolerance reproduces the small-signal resolver dead zone");

        SerialLinkParameters selectedParameters = legacyParameters;
        selectedParameters.resolvedRateQpStepSizeTolerance = 1e-6;
        selectedParameters.mpccQpStepSizeTolerance = 1e-4;
        check(selectedParameters.resolved_rate_qp_options().stepSizeTolerance == 1e-6
              && selectedParameters.mpcc_qp_options().stepSizeTolerance == 1e-4
              && selectedParameters.qpsolver.stepSizeTolerance == 5e-2,
              "split resolver/MPCC tolerances preserve the legacy fallback value");
        auto selectedModel = make_model(fixture.path);
        RobotLibrary::Control::SerialLinkKinematic selected(
            selectedModel, "tool", selectedParameters);
        const Eigen::VectorXd selectedCommand = selected.resolve_endpoint_twist(request);
        const Eigen::Vector<double,6> selectedTwist =
            selectedModel->jacobian(selectedModel->find_frame("tool")) * selectedCommand;
        check((selectedTwist - request).norm() / request.norm() <= 0.02,
              "selected resolver tolerance reproduces a 0.05 mm/s request within 2 percent");

        SerialLinkParameters referenceParameters = selectedParameters;
        referenceParameters.resolvedRateQpStepSizeTolerance = 1e-8;
        auto referenceModel = make_model(fixture.path);
        RobotLibrary::Control::SerialLinkKinematic reference(
            referenceModel, "tool", referenceParameters);
        const Eigen::VectorXd referenceCommand = reference.resolve_endpoint_twist(request);
        check((selectedCommand - referenceCommand).norm() <= 1e-10,
              "selected resolver solution matches the high-accuracy reference");

        SerialLinkParameters singularParameters = selectedParameters;
        singularParameters.minManipulability = 1.0;
        auto singularModel = make_model(fixture.path);
        RobotLibrary::Control::SerialLinkKinematic singular(
            singularModel, "tool", singularParameters);
        check(singular.resolve_endpoint_twist(request).allFinite(),
              "selected tolerance leaves the damped singular branch finite");
    }

    // The selected outer tolerance converges to the same physical first command
    // as 1e-8, and its frozen condensed inequalities remain satisfied.
    {
        SerialLinkParameters selectedParameters;
        selectedParameters.qpsolver.stepSizeTolerance = 5e-2;
        selectedParameters.mpccQpStepSizeTolerance = 1e-4;
        selectedParameters.qpsolver.maxSteps = 100;
        SerialLinkParameters referenceParameters = selectedParameters;
        referenceParameters.mpccQpStepSizeTolerance = 1e-8;
        auto selectedModel = make_model(fixture.path);
        auto referenceModel = make_model(fixture.path);
        SerialLinkMPCC selected(
            selectedModel, "tool", selectedParameters, 12, 0.002);
        SerialLinkMPCC reference(
            referenceModel, "tool", referenceParameters, 12, 0.002);
        set_moving_path(selected);
        set_moving_path(reference);
        (void)selected.step_at_time(0.0, 0.002);
        (void)reference.step_at_time(0.0, 0.002);
        const auto &selectedDiagnostics = selected.diagnostics();
        const auto &referenceDiagnostics = reference.diagnostics();
        check((selectedDiagnostics.optimalHorizon.head<3>()
               - referenceDiagnostics.optimalHorizon.head<3>()).norm() <= 1e-6,
              "selected MPCC tolerance converges within 1 um/s of reference");
        check((selectedDiagnostics.qpConstraintMatrix
               * selectedDiagnostics.optimalHorizon
               - selectedDiagnostics.qpConstraintVector).maxCoeff() <= 1e-8,
              "selected MPCC solution satisfies frozen inequalities");
    }

    if(failures == 0)
    {
        std::cout << "maximum affine-vs-explicit force error: "
                  << maximumAffineForceError << " N\n";
        std::cout << "predictive_contact_mpcc_test PASS\n";
        return 0;
    }
    std::cerr << "predictive_contact_mpcc_test FAIL (" << failures << " checks)\n";
    return 1;
}
