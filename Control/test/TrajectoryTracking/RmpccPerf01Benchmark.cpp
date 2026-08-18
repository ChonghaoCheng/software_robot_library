/**
 * @file RmpccPerf01Benchmark.cpp
 * @brief Deterministic offline solve_rmpcc component timing sweep for PERF-01.
 */

#include <Control/TrajectoryTracking/SerialLinkRMPCC.h>
#include <Model/KinematicTree.h>

#include <Eigen/Geometry>

#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <system_error>
#include <vector>

namespace {

using RobotLibrary::Control::RmpccParameters;
using RobotLibrary::Control::SerialLinkParameters;
using RobotLibrary::Control::SerialLinkRMPCC;
using RobotLibrary::Model::Pose;
using RobotLibrary::Trajectory::CartesianSpline;

class UrdfFixture
{
public:
    UrdfFixture()
    {
        const auto nonce = std::chrono::steady_clock::now()
                               .time_since_epoch().count();
        path = std::filesystem::temp_directory_path()
            / ("rmpcc_perf01_" + std::to_string(nonce) + ".urdf");
        std::ofstream stream(path);
        stream << R"(<robot name="perf01_fixture">
  <link name="base"/><link name="arm"/><link name="tool"/>
  <joint name="joint" type="revolute">
    <parent link="base"/><child link="arm"/><axis xyz="0 0 1"/>
    <limit lower="-3.14" upper="3.14" effort="10" velocity="1"/>
  </joint>
  <joint name="tool_fixed" type="fixed">
    <parent link="arm"/><child link="tool"/>
    <origin xyz="0.1 0 0" rpy="0 0 0"/>
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

CartesianSpline make_path()
{
    std::vector<Pose> poses;
    std::vector<double> times;
    for(int i = 0; i < 9; ++i)
    {
        const double u = static_cast<double>(i) / 8.0;
        const Eigen::Vector3d p(
            0.10 + 0.08 * u,
            0.025 * std::sin(2.0 * M_PI * u),
            0.020 * std::cos(M_PI * u));
        const Eigen::AngleAxisd rx(0.8 * u, Eigen::Vector3d::UnitX());
        const Eigen::AngleAxisd ry(0.5 * std::sin(M_PI * u), Eigen::Vector3d::UnitY());
        const Eigen::AngleAxisd rz(1.2 * u, Eigen::Vector3d::UnitZ());
        poses.emplace_back(p, Eigen::Quaterniond(rz * ry * rx));
        times.push_back(6.0 * u);
    }
    return CartesianSpline(poses, times, Eigen::Vector<double,6>::Zero());
}

} // namespace

int main(int argc, char **argv)
{
    int repetitions = 11;
    if(argc == 2)
    {
        repetitions = std::stoi(argv[1]);
    }
    if(repetitions < 1)
    {
        std::fprintf(stderr, "repetitions must be positive\n");
        return 2;
    }

    UrdfFixture fixture;
    auto model = std::make_shared<RobotLibrary::Model::KinematicTree>(
        fixture.path.string());
    model->update_state(Eigen::VectorXd::Zero(1), Eigen::VectorXd::Zero(1));
    const CartesianSpline path = make_path();
    constexpr double dt = 0.002;
    constexpr double progress = 0.30;

    std::cout << "horizon,repetition,reference_pose_seconds,reference_pose_count,reference_pose_requests,"
                 "reference_tangent_seconds,reference_tangent_count,reference_tangent_requests,"
                 "rollout_linearization_seconds,residual_hessian_seconds,"
                 "constraint_seconds,qp_seconds,total_seconds,"
                 "u0,u1,u2,u3,u4,u5,sdot,state_hash,residual_hash,hessian_hash\n";
    std::cout << std::setprecision(17);
    for(const int horizon : {20, 30, 40, 50, 75, 100})
    {
        RmpccParameters parameters;
        parameters.horizonSteps = horizon;
        SerialLinkRMPCC controller(
            model, "tool", SerialLinkParameters{}, parameters);
        controller.set_trajectory(path);
        for(int warmup = 0; warmup < 3; ++warmup)
        {
            controller.step(dt, progress);
        }
        for(int repetition = 0; repetition < repetitions; ++repetition)
        {
            controller.step(dt, progress);
            const auto &d = controller.diagnostics();
            std::cout << horizon << ',' << repetition << ','
                      << d.referencePoseQueryTimeSeconds << ','
                      << d.referencePoseQueryCount << ','
                      << d.referencePoseRequestCount << ','
                      << d.referenceTangentQueryTimeSeconds << ','
                      << d.referenceTangentQueryCount << ','
                      << d.referenceTangentRequestCount << ','
                      << d.stageRolloutLinearizationTimeSeconds << ','
                      << d.residualHessianAssemblyTimeSeconds << ','
                      << d.constraintConstructionTimeSeconds << ','
                      << d.qpSolveTimeSeconds << ','
                      << d.totalSolveTimeSeconds;
            for(int i = 0; i < 6; ++i)
            {
                std::cout << ',' << d.bodyTwist(i);
            }
            std::cout << ',' << d.progressRate << ','
                      << d.stateLinearizationHash << ','
                      << d.residualLinearizationHash << ','
                      << d.hessianHash << '\n';
        }
    }
    return 0;
}
