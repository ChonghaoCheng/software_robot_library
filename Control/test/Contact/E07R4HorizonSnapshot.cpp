/** E07-R4 non-production full predictive-MPCC QP snapshot utility. */

#include <Control/Contact/SerialLinkPredictiveContactMPCC.h>
#include <Model/KinematicTree.h>
#include <Trajectory/CartesianSpline.h>

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace {

using RobotLibrary::Control::MpccObjectiveParameters;
using RobotLibrary::Control::PredictiveContactMode;
using RobotLibrary::Control::PredictiveContactMpccParameters;
using RobotLibrary::Control::SerialLinkParameters;
using RobotLibrary::Control::SerialLinkPredictiveContactMPCC;
using RobotLibrary::Model::KinematicTree;
using RobotLibrary::Model::Pose;
using RobotLibrary::Trajectory::CartesianSpline;
using RobotLibrary::Trajectory::CartesianTrajectoryFrameState;

constexpr double kDt = 0.002;
constexpr double kDuration = 10.0;
constexpr double kAmplitude = 25e-6;
constexpr double kFrequency = 0.7;

void write_matrix(const std::filesystem::path &path,
                  const Eigen::MatrixXd &matrix)
{
    std::ofstream stream(path);
    stream << std::setprecision(17);
    for(Eigen::Index row = 0; row < matrix.rows(); ++row)
    {
        for(Eigen::Index column = 0; column < matrix.cols(); ++column)
        {
            if(column) stream << ',';
            stream << matrix(row, column);
        }
        stream << '\n';
    }
}

Eigen::Matrix4d board_pose(const bool moving, const double elapsed)
{
    Eigen::Matrix4d result = Eigen::Matrix4d::Identity();
    result.block<3,1>(0,3) = Eigen::Vector3d(-0.154, -0.454, 0.388);
    if(moving)
        result(1,3) += kAmplitude * std::sin(2.0 * M_PI * kFrequency * elapsed);
    return result;
}

double number(const char *text) { return std::stod(text); }

} // namespace

int main(int argc, char **argv)
{
    // urdf out horizon force z history[4] previous progress condition elapsed
    // board-frame position error xyz
    if(argc != 17)
    {
        std::cerr << "usage: e07_r4_horizon_snapshot URDF OUT N FORCE Z "
                     "H0 H1 H2 H3 PREVIOUS PROGRESS CONDITION ELAPSED EX EY EZ\n";
        return 2;
    }
    const std::string urdf = argv[1];
    const std::filesystem::path output = argv[2];
    const int horizon = std::stoi(argv[3]);
    const double force = number(argv[4]);
    const double realization = number(argv[5]);
    std::vector<double> history;
    for(int index = 6; index < 10; ++index) history.push_back(number(argv[index]));
    const double previous = number(argv[10]);
    const double progress = number(argv[11]);
    const bool moving = std::string(argv[12]) == "V1";
    const double elapsed = number(argv[13]);
    const Eigen::Vector3d loggedError(number(argv[14]), number(argv[15]),
                                      number(argv[16]));
    (void)previous; // The final history entry is the issued previous command.

    std::filesystem::create_directories(output);
    const std::array<double,6> qRaw{
        1.5708, -1.5708, -1.5708, -3.14, -1.5708, 0.0};
    const Eigen::VectorXd q = Eigen::Map<const Eigen::Vector<double,6>>(qRaw.data());
    auto model = std::make_shared<KinematicTree>(urdf);
    model->update_state(q, Eigen::VectorXd::Zero(6));

    SerialLinkParameters solver;
    solver.controlFrequency = 500.0;
    solver.minManipulability = 0.001;
    solver.maxJointAcceleration = 10.0;
    solver.qpsolver.method = "active set";
    solver.qpsolver.maxSteps = 50;
    solver.resolvedRateQpStepSizeTolerance = 1e-6;
    solver.mpccQpStepSizeTolerance = 1e-4;
    SerialLinkPredictiveContactMPCC controller(
        model, "wrist_3_link", solver, static_cast<unsigned int>(horizon), kDt);

    PredictiveContactMpccParameters contact;
    contact.mode = PredictiveContactMode::Active;
    contact.compressionDirectionParent = Eigen::Vector3d(0.0, -1.0, 0.0);
    contact.contactOffsetEndpoint = Eigen::Vector3d(
        0.0005265895, 0.0999391081, -0.0399956612);
    contact.stiffness = 21303.75539503847;
    contact.targetForce = 2.5;
    contact.forceWeight = 1.0;
    contact.minimumForce = 1.8;
    contact.maximumForce = 3.2;
    contact.maximumPenetrationIncrement = 4e-6;
    contact.forceSlackWeight = 1e5;
    contact.penetrationSlackWeight = 1e16;
    contact.actuationAware = true;
    contact.realizationAutoregressive = 0.9809437967444158;
    contact.realizationInputGain = 0.007248865304869999;
    contact.realizationDelay = 4;
    contact.minimumContactModelForce = 0.5;
    contact.normalActionGuardEnabled = true;
    contact.maximumRobotNormalCommand = 1e-4;
    contact.maximumRobotNormalCommandStep = 1e-6;
    contact.normalCommandSmoothWeight = 1e-6;
    controller.set_predictive_contact_parameters(contact);
    controller.set_linear_velocity_limit(0.18);
    controller.set_angular_velocity_limit(1.20);

    MpccObjectiveParameters objective;
    objective.contourWeight = 260.0;
    objective.lagWeight = 35.0;
    objective.orientationWeight = 40.0;
    objective.inputLinearWeight = 0.02;
    objective.inputAngularWeight = 0.02;
    objective.inputProgressWeight = 2.0;
    objective.inputDifferenceWeight = 3e-5;
    objective.pathVelocityWeight = 0.0;
    // Runtime normalization: progressWeight * nominalRate * controllerRate.
    objective.progressReward = 2.0 * 0.1 * 500.0;
    controller.set_objective_parameters(objective);
    controller.set_progress_rate_limits(0.0, 0.1, 0.125);

    const Eigen::Matrix4d currentBoard = board_pose(moving, elapsed);
    const Pose endpoint = controller.endpoint_pose();
    const Eigen::Matrix3d boardRotation = currentBoard.block<3,3>(0,0);
    const Eigen::Vector3d endpointBoard = boardRotation.transpose()
        * (endpoint.translation() - currentBoard.block<3,1>(0,3));
    const Eigen::Vector3d referenceAtProgress = endpointBoard - loggedError;
    const Eigen::Vector3d pathStart = referenceAtProgress
        - Eigen::Vector3d(0.03 * progress, 0.0, 0.0);
    const Eigen::Quaterniond orientationBoard(
        boardRotation.transpose() * endpoint.quaternion().toRotationMatrix());
    const Pose start(pathStart, orientationBoard);
    const Pose end(pathStart + Eigen::Vector3d(0.03, 0.0, 0.0), orientationBoard);
    controller.set_trajectory(CartesianSpline(
        start, end, Eigen::Vector<double,6>::Zero(), 0.0, kDuration));

    CartesianTrajectoryFrameState frame;
    if(moving)
    {
        frame.transformInBase = board_pose(true, elapsed - kDt);
        controller.set_trajectory_frame(frame, elapsed - kDt);
    }
    frame.transformInBase = currentBoard;
    controller.set_trajectory_frame(frame, elapsed);
    controller.set_force_measurement(force, true);
    controller.set_normal_realization_state(realization, history);

    std::string exception;
    try { (void)controller.step_at_time(elapsed, kDt, progress); }
    catch(const std::exception &error) { exception = error.what(); }
    const auto &diagnostics = controller.diagnostics();
    const auto &contactDiagnostics = controller.predictive_contact_diagnostics();
    write_matrix(output / "H.csv", diagnostics.qpHessian);
    write_matrix(output / "f.csv", diagnostics.qpGradient);
    write_matrix(output / "B.csv", diagnostics.qpConstraintMatrix);
    write_matrix(output / "z.csv", diagnostics.qpConstraintVector);
    // Fixed variables are carried as an equality block, so a snapshot that
    // records only B/z is no longer the complete frozen QP.
    write_matrix(output / "Aeq.csv", diagnostics.qpEqualityMatrix);
    write_matrix(output / "yeq.csv", diagnostics.qpEqualityVector);
    write_matrix(output / "U0.csv", diagnostics.qpSeed);
    write_matrix(output / "production_solution.csv",
                 diagnostics.qpReturnedSolution);
    write_matrix(output / "predicted_force.csv",
                 contactDiagnostics.predictedForce);
    write_matrix(output / "predicted_robot_normal_command.csv",
                 contactDiagnostics.predictedCommandedRobotNormalVelocity);
    std::ofstream metadata(output / "metadata.json");
    metadata << std::setprecision(17)
             << "{\n  \"horizon\": " << horizon
             << ",\n  \"force_n\": " << force
             << ",\n  \"realization_z_mps\": " << realization
             << ",\n  \"progress\": " << progress
             << ",\n  \"condition\": \"" << (moving ? "V1" : "V0") << "\""
             << ",\n  \"elapsed_s\": " << elapsed
             << ",\n  \"production_qp_iterations\": "
             << diagnostics.qpIterations
             << ",\n  \"production_qp_solve_time_ms\": "
             << 1000.0 * diagnostics.qpSolveTimeSeconds
             << ",\n  \"production_qp_objective\": "
             << diagnostics.qpObjective
             << ",\n  \"production_qp_max_violation\": "
             << diagnostics.qpMaximumConstraintViolation
             << ",\n  \"force_cost\": " << contactDiagnostics.forceCost
             << ",\n  \"force_slack_max_n\": "
             << contactDiagnostics.maximumForceSlack
             << ",\n  \"penetration_slack_max_m\": "
             << contactDiagnostics.maximumPenetrationSlack
             << ",\n  \"runtime_exception\": \"" << exception << "\"\n}\n";
    std::cout << "wrote " << diagnostics.qpHessian.rows() << " variables, "
              << diagnostics.qpConstraintMatrix.rows() << " constraints\n";
    return diagnostics.qpHessian.size() == 0 ? 1 : 0;
}
