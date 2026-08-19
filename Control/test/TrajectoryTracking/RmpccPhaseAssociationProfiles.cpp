/**
 * @file TrajectoryTracking/test/RmpccPhaseAssociationProfiles.cpp
 * @brief Permanent constructor gates for phase-association/profile compatibility.
 */

#include <Control/TrajectoryTracking/SerialLinkRMPCC.h>
#include <Control/TrajectoryTracking/SerialLinkMPCC.h>
#include <Model/KinematicTree.h>

#include <Eigen/Core>

#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <system_error>

namespace {

using RobotLibrary::Control::RmpccContourResidualGeometry;
using RobotLibrary::Control::RmpccLagGeometry;
using RobotLibrary::Control::RmpccLagPenalty;
using RobotLibrary::Control::RmpccObjectiveGeometry;
using RobotLibrary::Control::RmpccParameters;
using RobotLibrary::Control::RmpccPhaseAssociation;
using RobotLibrary::Control::SerialLinkParameters;
using RobotLibrary::Control::SerialLinkRMPCC;

class UrdfFixture
{
public:
    UrdfFixture()
    {
        const auto nonce = std::chrono::steady_clock::now()
                               .time_since_epoch().count();
        path = std::filesystem::temp_directory_path()
            / ("rmpcc_phase_profiles_" + std::to_string(nonce) + ".urdf");
        std::ofstream stream(path);
        stream << R"(<robot name="phase_profile_fixture">
  <link name="base"/>
  <link name="arm"/>
  <link name="tool"/>
  <joint name="joint" type="revolute">
    <parent link="base"/>
    <child link="arm"/>
    <axis xyz="0 0 1"/>
    <limit lower="-3.14" upper="3.14" effort="10" velocity="1"/>
  </joint>
  <joint name="tool_fixed" type="fixed">
    <parent link="arm"/>
    <child link="tool"/>
    <origin xyz="0.1 0 0" rpy="0 0 0"/>
  </joint>
</robot>)";
        stream.close();
    }

    ~UrdfFixture()
    {
        std::error_code error;
        std::filesystem::remove(path, error);
    }

    std::filesystem::path path;
};

RmpccParameters local_unified(const RmpccPhaseAssociation association)
{
    RmpccParameters value;
    value.phaseAssociation = association;
    value.rotationCharacteristicLength = 0.20;
    return value;
}

RmpccParameters associated(
    const RmpccPhaseAssociation association,
    const RmpccContourResidualGeometry geometry)
{
    RmpccParameters value = local_unified(association);
    value.contourResidualGeometry = geometry;
    value.lagPenalty = RmpccLagPenalty::ScalarPosePathArc;
    if(geometry
       == RmpccContourResidualGeometry::AssociatedDecoupledCartesianSO3)
    {
        value.objectiveGeometry = RmpccObjectiveGeometry::DecoupledCartesianSO3;
    }
    return value;
}

RmpccParameters inert_scheduled_decoupled()
{
    RmpccParameters value = local_unified(
        RmpccPhaseAssociation::TaskPoseFeature);
    value.objectiveGeometry = RmpccObjectiveGeometry::DecoupledCartesianSO3;
    value.contourResidualGeometry =
        RmpccContourResidualGeometry::ScheduledDecoupledCartesianSO3;
    value.runningLagGeometry = RmpccLagGeometry::SplitTranslationRotation;
    return value;
}

bool accepts(const std::shared_ptr<RobotLibrary::Model::KinematicTree> &model,
             const RmpccParameters &profile,
             const char *name)
{
    try
    {
        SerialLinkRMPCC controller(
            model, "tool", SerialLinkParameters{}, profile);
        return true;
    }
    catch(const std::exception &error)
    {
        std::fprintf(stderr, "%s unexpectedly rejected: %s\n", name, error.what());
        return false;
    }
}

bool rejects(const std::shared_ptr<RobotLibrary::Model::KinematicTree> &model,
             const RmpccParameters &profile,
             const char *expectedMessage,
             const char *name)
{
    try
    {
        SerialLinkRMPCC controller(
            model, "tool", SerialLinkParameters{}, profile);
    }
    catch(const std::invalid_argument &error)
    {
        if(std::string(error.what()).find(expectedMessage) != std::string::npos)
        {
            return true;
        }
        std::fprintf(stderr, "%s rejected for the wrong reason: %s\n",
                     name, error.what());
        return false;
    }
    catch(const std::exception &error)
    {
        std::fprintf(stderr, "%s threw the wrong exception: %s\n",
                     name, error.what());
        return false;
    }
    std::fprintf(stderr, "%s unexpectedly accepted\n", name);
    return false;
}

} // namespace

int main()
{
    UrdfFixture fixture;
    auto model = std::make_shared<RobotLibrary::Model::KinematicTree>(
        fixture.path.string());
    model->update_state(Eigen::VectorXd::Zero(1), Eigen::VectorXd::Zero(1));

    int failures = 0;
    {
        RobotLibrary::Control::SerialLinkMPCC controller(
            model, "tool", SerialLinkParameters{});
        failures += controller.linear_velocity_limit() != 0.2;
        failures += controller.angular_velocity_limit() != 0.2;
        controller.set_angular_velocity_limit(0.5);
        failures += controller.angular_velocity_limit() != 0.5;
        for(const double invalidLimit :
            {0.0, -0.5, std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::quiet_NaN()})
        {
            try
            {
                controller.set_angular_velocity_limit(invalidLimit);
                ++failures;
            }
            catch(const std::invalid_argument &)
            {
            }
        }
    }
    failures += not accepts(
        model, local_unified(RmpccPhaseAssociation::TaskPoseFeature),
        "TaskPoseFeature + LocalUnified");
    failures += not rejects(
        model, inert_scheduled_decoupled(), "has no cost effect",
        "TaskPoseFeature + inert scheduled decoupled");
    for(const auto geometry :
        {RmpccContourResidualGeometry::ExactAssociatedUnifiedSE3,
         RmpccContourResidualGeometry::AssociatedDecoupledCartesianSO3})
    {
        failures += not accepts(
            model, associated(RmpccPhaseAssociation::TaskPointXYZ, geometry),
            "associated + TaskPointXYZ");
        failures += not accepts(
            model, associated(RmpccPhaseAssociation::TaskPoseFeature, geometry),
            "associated + TaskPoseFeature");
        failures += not rejects(
            model, associated(RmpccPhaseAssociation::MetricScrew, geometry),
            "require a task phase association", "associated + MetricScrew");
    }

    for(const double invalidLength :
        {0.0, -0.20, std::numeric_limits<double>::infinity(),
         std::numeric_limits<double>::quiet_NaN()})
    {
        RmpccParameters profile = local_unified(
            RmpccPhaseAssociation::TaskPoseFeature);
        profile.rotationCharacteristicLength = invalidLength;
        failures += not rejects(
            model, profile, "requires a positive rotationCharacteristicLength",
            "TaskPoseFeature without positive l_R");
    }

    return failures == 0 ? 0 : 1;
}
