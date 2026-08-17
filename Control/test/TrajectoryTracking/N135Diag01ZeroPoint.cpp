/**
 * @file N135Diag01ZeroPoint.cpp
 * @brief DIAG-01 exact-on-path path-velocity identity gate.
 */

#include "RmpccAssociatedPhase.h"
#include "RmpccCostGeometry.h"
#include "RmpccPhaseResidual.h"
#include "RmpccPrediction.h"
#include "RmpccResidualLinearization.h"

#include <Math/MathFunctions.h>
#include <Model/Pose.h>
#include <Trajectory/CartesianSpline.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using RobotLibrary::Control::RmpccAssociatedResiduals;
using RobotLibrary::Control::RmpccContourResidualGeometry;
using RobotLibrary::Control::RmpccPhaseAssociation;
using RobotLibrary::Control::RmpccPoseArcTable;
using RobotLibrary::Control::RmpccStateVector;
using RobotLibrary::Model::Pose;
using RobotLibrary::Trajectory::CartesianSpline;

struct Path
{
    std::vector<Pose> poses;
    std::vector<double> times;
};

struct Profile
{
    const char *name;
    const char *phase;
    const char *objective;
    bool associated;
    bool decoupled;
};

struct Result
{
    std::string frame;
    std::string profile;
    double progress = 0.0;
    double linear = 0.0;
    double angular = 0.0;
    double total = 0.0;
};

std::vector<std::string> split(const std::string &line)
{
    std::vector<std::string> result;
    std::stringstream stream(line);
    std::string item;
    while(std::getline(stream, item, ',')) result.push_back(item);
    return result;
}

Path load_path(const std::filesystem::path &file)
{
    std::ifstream stream(file);
    if(not stream) throw std::runtime_error("cannot open " + file.string());
    std::string line;
    std::getline(stream, line);
    const std::vector<std::string> header = split(line);
    const auto column = [&](const std::string &name)
    {
        for(std::size_t i = 0; i < header.size(); ++i)
            if(header[i] == name) return i;
        throw std::runtime_error("missing CSV column " + name);
    };
    const std::size_t time = column("time");
    const std::size_t x = column("x"), y = column("y"), z = column("z");
    const std::size_t roll = column("roll"), pitch = column("pitch");
    const std::size_t yaw = column("yaw");
    Path result;
    while(std::getline(stream, line))
    {
        if(line.empty()) continue;
        const std::vector<std::string> value = split(line);
        const Eigen::Vector3d position(
            std::stod(value[x]), std::stod(value[y]), std::stod(value[z]));
        const Eigen::Quaterniond orientation(
            Eigen::AngleAxisd(std::stod(value[yaw]), Eigen::Vector3d::UnitZ())
            * Eigen::AngleAxisd(std::stod(value[pitch]), Eigen::Vector3d::UnitY())
            * Eigen::AngleAxisd(std::stod(value[roll]), Eigen::Vector3d::UnitX()));
        result.poses.emplace_back(position, orientation);
        result.times.push_back(std::stod(value[time]));
    }
    if(result.poses.size() < 2) throw std::runtime_error("path has fewer than two poses");
    return result;
}

Eigen::Matrix<double,6,6> lag_weight()
{
    Eigen::Matrix<double,6,6> value = Eigen::Matrix<double,6,6>::Zero();
    value.diagonal() << 10.0, 10.0, 10.0, 0.4, 0.4, 0.4;
    return value;
}

template<typename ReferenceTransformFunction, typename ReferenceTangentFunction>
void exercise_branch(const Profile &profile, const RmpccStateVector &state,
                     const RmpccPoseArcTable &arc,
                     ReferenceTransformFunction &&transform,
                     ReferenceTangentFunction &&tangent)
{
    if(profile.associated)
    {
        const auto geometry = profile.decoupled
            ? RmpccContourResidualGeometry::AssociatedDecoupledCartesianSO3
            : RmpccContourResidualGeometry::ExactAssociatedUnifiedSE3;
        const RmpccAssociatedResiduals residual =
            RobotLibrary::Control::rmpcc_associated_residuals(
                state, geometry, 1e-12, arc, transform, tangent);
        if(not residual.context.observable)
            throw std::runtime_error(std::string(profile.name) + " is unobservable");
    }
    else if(profile.decoupled)
    {
        const auto residual = RobotLibrary::Control::rmpcc_decoupled_residuals(
            state, 1e-8, tangent);
        if(not residual.contour.allFinite())
            throw std::runtime_error(std::string(profile.name) + " is non-finite");
    }
    else
    {
        Eigen::Matrix<double,6,6> metric = Eigen::Matrix<double,6,6>::Identity();
        metric.diagonal().tail<3>().setConstant(0.04);
        const auto association = std::string(profile.phase) == "task"
            ? RmpccPhaseAssociation::TaskPointXYZ
            : RmpccPhaseAssociation::MetricScrew;
        const auto residual = RobotLibrary::Control::rmpcc_phase_residuals(
            state, metric, association, 1e-8, 1e-12, transform, tangent);
        if(not residual.phaseObservable)
            throw std::runtime_error(std::string(profile.name) + " is unobservable");
    }
}

} // namespace

int main(int argc, char **argv)
{
    try
    {
        if(argc != 4)
            throw std::invalid_argument(
                "usage: n135_diag01_zero_point D000_CSV D050_CSV OUTPUT_YAML");
        const std::array<std::pair<const char *, std::filesystem::path>,2> paths{{
            {"d000", argv[1]}, {"d050", argv[2]}}};
        const std::array<Profile,6> profiles{{
            {"MV", "metric", "unified", false, false},
            {"TV", "task", "unified", false, false},
            {"TG", "task", "decoupled_vector", true, false},
            {"TPA", "task", "decoupled_arc", true, false},
            {"DS", "scheduled", "decoupled", false, true},
            {"DT", "task", "decoupled", true, true}}};
        const std::array<double,10> progress{{
            0.05, 0.15, 0.25, 0.35, 0.45,
            0.55, 0.65, 0.75, 0.85, 0.95}};
        constexpr double timeScale = 6.0;
        constexpr double tolerance = 1e-10;
        std::vector<Result> results;
        double scheduleRate = 0.0;
        double maximum = 0.0;

        for(const auto &[frame, file] : paths)
        {
            Path input = load_path(file);
            for(double &time : input.times) time *= timeScale;
            CartesianSpline trajectory(
                input.poses, input.times, Eigen::Vector<double,6>::Zero());
            scheduleRate = 1.0 / (trajectory.end_time() - trajectory.start_time());
            const auto transform = [&](const double s)
            {
                return trajectory.pose_at_progress(s).as_matrix();
            };
            const auto tangent = [&](const double s)
            {
                return trajectory.tangent_at_progress(s, 1e-3);
            };
            RmpccPoseArcTable arc;
            arc.build(lag_weight(), tangent, 4097);
            for(const Profile &profile : profiles)
            {
                for(const double s : progress)
                {
                    const Eigen::Matrix4d reference = transform(s);
                    const Eigen::Matrix4d actual = reference;
                    RmpccStateVector state = RmpccStateVector::Zero();
                    state(6) = s;
                    state.head<6>() = RobotLibrary::Math::se3_logarithm(
                        RobotLibrary::Math::se3_inverse(reference) * actual);
                    exercise_branch(profile, state, arc, transform, tangent);
                    const Eigen::Vector<double,6> tau = tangent(s);
                    const Eigen::Vector<double,6> command = tau * scheduleRate;
                    const Eigen::Vector<double,6> residual = command
                        - RobotLibrary::Control::rmpcc_transport_reference_tangent(
                            state.head<6>(), tau) * scheduleRate;
                    Result result{frame, profile.name, s,
                                  residual.head<3>().norm(),
                                  residual.tail<3>().norm(), residual.norm()};
                    maximum = std::max(maximum, result.total);
                    results.push_back(result);
                    if(result.total >= tolerance)
                        throw std::runtime_error(
                            result.profile + " zero-point identity failed");
                }
            }
        }

        std::ofstream output(argv[3]);
        if(not output) throw std::runtime_error("cannot open output YAML");
        output << std::setprecision(17)
               << "schema: n135_diag01_zero_point_v1\n"
               << "experiment_id: N135\n"
               << "stage: A1\n"
               << "trajectory_family: board_se3_bridge\n"
               << "task_frames: [d000, d050]\n"
               << "trajectory_time_scale: " << timeScale << "\n"
               << "schedule_rate: " << scheduleRate << "\n"
               << "tolerance: " << tolerance << "\n"
               << "maximum_residual: " << maximum << "\n"
               << "all_passed: true\n"
               << "a2_a3_skipped: true\n"
               << "h1: falsified\n"
               << "samples:\n";
        for(const Result &result : results)
            output << "  - {frame: " << result.frame
                   << ", profile: " << result.profile
                   << ", s: " << result.progress
                   << ", linear: " << result.linear
                   << ", angular: " << result.angular
                   << ", total: " << result.total << "}\n";
        return 0;
    }
    catch(const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
#include <algorithm>
