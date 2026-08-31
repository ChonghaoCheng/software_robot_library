#include <Control/TrajectoryTracking/ParentFrameReferenceMotion.h>
#include <Control/TrajectoryTracking/SerialLinkLieAlgebraMPC.h>
#include <Control/TrajectoryTracking/SerialLinkMPCC.h>
#include <Control/TrajectoryTracking/SerialLinkRMPCC.h>
#include <Math/MathFunctions.h>
#include <Model/KinematicTree.h>
#include <Trajectory/CartesianSpline.h>

#include <Eigen/Geometry>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <thread>
#include <chrono>

namespace {

int failures = 0;

void check(bool condition, const std::string &message)
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
            / ("causal_geometry_prediction_" + std::to_string(nonce) + ".urdf");
        std::ofstream stream(path);
        stream << R"(<robot name="prediction_fixture">
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
    ~UrdfFixture() { std::error_code error; std::filesystem::remove(path, error); }
    std::filesystem::path path;
};

RobotLibrary::Trajectory::CartesianSpline make_path(
    const RobotLibrary::Model::Pose &start)
{
    std::vector<RobotLibrary::Model::Pose> poses;
    std::vector<double> times;
    for(int i = 0; i < 9; ++i)
    {
        const double u = static_cast<double>(i) / 8.0;
        poses.emplace_back(
            start.translation() + Eigen::Vector3d(0.03*u, 0.01*std::sin(M_PI*u), 0.01*u),
            Eigen::Quaterniond(Eigen::AngleAxisd(0.2*u, Eigen::Vector3d::UnitY()))
                * start.quaternion());
        times.push_back(4.0*u);
    }
    return RobotLibrary::Trajectory::CartesianSpline(
        poses, times, Eigen::Vector<double,6>::Zero());
}

std::shared_ptr<RobotLibrary::Model::KinematicTree> model_at_state(
    const std::filesystem::path &path)
{
    auto model = std::make_shared<RobotLibrary::Model::KinematicTree>(path.string());
    model->update_state(
        Eigen::VectorXd::Constant(model->number_of_joints(), 0.25),
        Eigen::VectorXd::Zero(model->number_of_joints()));
    return model;
}

void prediction_kinematics()
{
    using RobotLibrary::Control::CausalParentFrameMotion;
    using RobotLibrary::Control::predicted_parent_frame_state;
    using RobotLibrary::Math::se3_exponential;
    using RobotLibrary::Math::se3_inverse;
    using RobotLibrary::Math::se3_logarithm;

    const double sampleDt = 0.01;
    const double horizonDt = 0.002;
    const Eigen::Matrix4d initial =
        se3_exponential((Eigen::Vector<double,6>() <<
            0.03, -0.02, 0.04, 0.1, -0.07, 0.05).finished());
    const std::vector<Eigen::Vector<double,6>> cases = {
        (Eigen::Vector<double,6>() << 0.08, -0.03, 0.04, 0, 0, 0).finished(),
        (Eigen::Vector<double,6>() << 0, 0, 0, 0.12, -0.09, 0.07).finished(),
        (Eigen::Vector<double,6>() << 0.05, 0.02, -0.01, -0.08, 0.1, 0.04).finished()};

    for(std::size_t c = 0; c < cases.size(); ++c)
    {
        CausalParentFrameMotion motion;
        motion.update(initial, 1.0);
        const Eigen::Matrix4d measuredNext = initial * se3_exponential(sampleDt * cases[c]);
        motion.update(measuredNext, 1.0 + sampleDt);
        check((motion.body_twist() - cases[c]).norm() < 2e-12,
              "two-sample body-twist estimate case " + std::to_string(c));
        for(int stage : {0, 1, 7, 20})
        {
            const Eigen::Matrix4d expected =
                measuredNext * se3_exponential(stage * horizonDt * cases[c]);
            check(se3_logarithm(se3_inverse(expected)
                                * motion.predicted_pose(stage, horizonDt)).norm() < 2e-12,
                  "predicted pose case " + std::to_string(c));
            const auto frame = predicted_parent_frame_state(motion, stage, horizonDt);
            const Eigen::Matrix3d rotation = expected.block<3,3>(0,0);
            check((frame.twistInBase.head<3>() - rotation * cases[c].head<3>()).norm() < 2e-12,
                  "predicted base linear point twist case " + std::to_string(c));
            check((frame.twistInBase.tail<3>() - rotation * cases[c].tail<3>()).norm() < 2e-12,
                  "predicted base angular point twist case " + std::to_string(c));
        }
    }

    // Host pacing is outside the estimator contract: identical pose/timestamp
    // sequences must produce bitwise-equivalent estimator outputs.
    CausalParentFrameMotion steady;
    CausalParentFrameMotion jittered;
    for(int index = 0; index < 8; ++index)
    {
        const double time = 2.0 + sampleDt * index;
        const Eigen::Matrix4d pose = initial * se3_exponential(
            sampleDt * index * cases.back());
        steady.update(pose, time);
        if(index % 2 == 0)
            std::this_thread::sleep_for(std::chrono::microseconds(100));
        jittered.update(pose, time);
        check(steady.current_pose().isApprox(jittered.current_pose(), 0.0)
              && steady.body_twist().isApprox(jittered.body_twist(), 0.0)
              && steady.current_time() == jittered.current_time(),
              "parent estimator host-pacing invariance");
    }
    const auto timestampedBefore = steady.timestamped_setter_count();
    const auto duplicateBefore = steady.duplicate_timestamp_count();
    steady.update(steady.current_pose(), steady.current_time());
    check(steady.timestamped_setter_count() == timestampedBefore + 1
          && steady.duplicate_timestamp_count() == duplicateBefore + 1,
          "duplicate timestamp observability");
    const auto outOfOrderBefore = steady.out_of_order_timestamp_count();
    steady.update(steady.current_pose(), steady.current_time() - sampleDt);
    check(steady.out_of_order_timestamp_count() == outOfOrderBefore + 1,
          "out-of-order timestamp observability");
    check(steady.too_small_interval_observable_without_policy_change(),
          "too-small interval is explicitly observable under the temporal policy");
}

void static_parent_controller_equivalence(const std::filesystem::path &urdf)
{
    using RobotLibrary::Control::RmpccParameters;
    using RobotLibrary::Control::SerialLinkLieAlgebraMPC;
    using RobotLibrary::Control::SerialLinkMPCC;
    using RobotLibrary::Control::SerialLinkParameters;
    using RobotLibrary::Control::SerialLinkRMPCC;
    using RobotLibrary::Trajectory::CartesianTrajectoryFrameState;

    SerialLinkParameters base;
    base.qpsolver.maxSteps = 50;
    base.qpsolver.stepSizeTolerance = 1e-5;
    CartesianTrajectoryFrameState frame;
    frame.transformInBase = Eigen::Matrix4d::Identity();

    auto lieStaticModel = model_at_state(urdf);
    auto liePredModel = model_at_state(urdf);
    SerialLinkLieAlgebraMPC lieStatic(lieStaticModel, "tool", base, 20, 0.002);
    SerialLinkLieAlgebraMPC liePred(liePredModel, "tool", base, 20, 0.002);
    lieStatic.set_feedback_bandwidth(2.0, 2.0);
    liePred.set_feedback_bandwidth(2.0, 2.0);
    const auto liePath = make_path(lieStatic.endpoint_pose());
    lieStatic.set_trajectory(liePath);
    liePred.set_trajectory(liePath);
    lieStatic.set_trajectory_frame(frame);
    liePred.set_trajectory_frame(frame, 0.0);
    liePred.set_trajectory_frame(frame, 0.01);
    const Eigen::VectorXd lieStaticCommand = lieStatic.track_endpoint_trajectory_at_time(0.25);
    const Eigen::VectorXd liePredCommand = liePred.track_endpoint_trajectory_at_time(0.25);
    check((lieStaticCommand - liePredCommand).norm() < 1e-10,
          "Lie static parent first command equivalence");
    check(liePred.diagnostics().parentFrameBodyTwist.norm() < 1e-12,
          "Lie static parent predicted twist is zero");

    auto mpccStaticModel = model_at_state(urdf);
    auto mpccPredModel = model_at_state(urdf);
    SerialLinkMPCC mpccStatic(mpccStaticModel, "tool", base);
    SerialLinkMPCC mpccPred(mpccPredModel, "tool", base);
    // This test qualifies parent-frame equivalence, not the optimized-progress
    // active-set path.  Use the deterministic fixed-progress path so a
    // degenerate zero-error task QP cannot mask the parent semantic check.
    mpccStatic.set_fixed_progress_schedule(true);
    mpccPred.set_fixed_progress_schedule(true);
    const auto mpccPath = make_path(mpccStatic.endpoint_pose());
    mpccStatic.set_trajectory(mpccPath);
    mpccPred.set_trajectory(mpccPath);
    mpccStatic.set_trajectory_frame(frame);
    mpccPred.set_trajectory_frame(frame, 0.0);
    mpccPred.set_trajectory_frame(frame, 0.01);
    const Eigen::VectorXd mpccStaticCommand = mpccStatic.step(0.002);
    const Eigen::VectorXd mpccPredCommand = mpccPred.step(0.002);
    check((mpccStaticCommand - mpccPredCommand).norm() < 1e-9,
          "MPCC static parent first command equivalence");
    check(mpccPred.diagnostics().parentFrameBodyTwist.norm() < 1e-12,
          "MPCC static parent predicted twist is zero");

    RmpccParameters rp;
    rp.horizonSteps = 20;
    rp.referenceMotion = RobotLibrary::Control::RmpccReferenceMotion::FiniteStageExact;
    rp.predictorGeometry = RobotLibrary::Control::RmpccPredictorGeometry::ExactSE3;
    rp.residualLinearization =
        RobotLibrary::Control::RmpccResidualLinearization::FullResidualJacobian;
    rp.progressRateMinMultiplier = 0.25;
    rp.progressRateMax = 0.0;
    rp.progressRateMaxMultiplier = 1.0;
    rp.progressReward = 0.2;
    rp.pathVelocityWeight += rp.controlWeight;
    rp.controlWeight.setZero();
    rp.enableDetailedDiagnostics = false;
    auto gpModel = model_at_state(urdf);
    auto noPredModel = model_at_state(urdf);
    SerialLinkRMPCC gp(gpModel, "tool", base, rp);
    SerialLinkRMPCC noPred(noPredModel, "tool", base, rp);
    gp.set_fixed_progress_schedule(true);
    noPred.set_fixed_progress_schedule(true);
    const auto rmpccPath = make_path(gp.endpoint_pose());
    gp.set_trajectory(rmpccPath);
    noPred.set_trajectory(rmpccPath);
    gp.set_trajectory_frame(frame, 0.0);
    gp.set_trajectory_frame(frame, 0.01);
    noPred.set_trajectory_frame(frame);
    const Eigen::VectorXd gpCommand = gp.step(0.002);
    const Eigen::VectorXd noPredCommand = noPred.step(0.002);
    check((gpCommand - noPredCommand).norm() < 1e-9,
          "GP versus no-parent-prediction static first command equivalence");
    check((gp.diagnostics().bodyTwist - noPred.diagnostics().bodyTwist).norm() < 1e-9,
          "GP versus no-parent-prediction static task command equivalence");
}

} // namespace

int main()
{
    prediction_kinematics();
    UrdfFixture fixture;
    static_parent_controller_equivalence(fixture.path);
    if(failures == 0)
    {
        std::cout << "causal_geometry_prediction_test PASS\n";
        return 0;
    }
    std::cerr << "causal_geometry_prediction_test FAIL (" << failures << ")\n";
    return 1;
}
