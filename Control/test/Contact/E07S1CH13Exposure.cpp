/** Offline S1C fixed-progress convergence exposure utility. */

#include <Control/Contact/SerialLinkContactMPCC.h>
#include <Model/KinematicTree.h>
#include <Trajectory/CartesianSpline.h>

#include <Eigen/Dense>

#include <array>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

namespace {

using namespace RobotLibrary::Control;
using RobotLibrary::Model::KinematicTree;
using RobotLibrary::Model::Pose;
using RobotLibrary::Trajectory::CartesianSpline;

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

} // namespace

int main(int argc, char **argv)
{
    if(argc != 6)
    {
        std::cerr << "usage: e07_s1c_h13_exposure URDF OUTPUT HORIZON MODE FIXED\n";
        return 2;
    }
    const std::string urdf = argv[1];
    const std::filesystem::path output = argv[2];
    const unsigned int horizon = static_cast<unsigned int>(std::stoul(argv[3]));
    const std::string modeName = argv[4];
    const bool fixed = std::string(argv[5]) == "fixed";
    const ContactMode mode = modeName == "Loss"
        ? ContactMode::Loss : ContactMode::Constraint;

    std::filesystem::create_directories(output);
    const std::array<double,6> qRaw{
        1.5708, -1.5708, -1.5708, -3.14, -1.5708, 0.0};
    const Eigen::VectorXd q =
        Eigen::Map<const Eigen::Vector<double,6>>(qRaw.data());
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

    SerialLinkContactMPCC controller(
        model, "wrist_3_link", solver, horizon, 0.002);
    controller.set_fixed_progress_schedule(fixed);
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
    objective.progressReward = 2.0 * 0.1 * 500.0;
    controller.set_objective_parameters(objective);

    ContactParameters contact;
    contact.mode = mode;
    contact.normalAxisInBoard = Eigen::Vector3d(0.0, -1.0, 0.0);
    contact.tangentAxisInBoard = Eigen::Vector3d::UnitX();
    contact.targetForce = 2.5;
    contact.forceResponseGain = 21303.75539503847;
    contact.forceTolerance = 0.7;
    contact.forceWeight = 1.0;
    contact.slackWeight = 1e5;
    contact.maxForceSlack = 50.0;
    contact.tangentPositionWeight = 260.0;
    contact.pathLagPositionWeight = 35.0;
    contact.normalPositionWeight = 0.0;
    contact.orientationWeight = 40.0;
    controller.set_contact_parameters(contact);

    ContactState state;
    state.inContact = true;
    state.normalBase = Eigen::Vector3d(0.0, -1.0, 0.0);
    state.measuredNormalForce = 2.5;
    controller.set_contact_state(state);

    const Pose start = controller.endpoint_pose();
    const Pose end(start.translation() + Eigen::Vector3d(0.03, 0.0, 0.0),
                   start.quaternion());
    controller.set_trajectory(CartesianSpline(
        start, end, Eigen::Vector<double,6>::Zero(), 0.0, 10.0));
    // set_trajectory derives nominal limits; apply the frozen E07 values after it,
    // matching the action-server construction order.
    controller.set_progress_rate_limits(0.0, 0.1, 0.125);

    std::string exception;
    try { (void)controller.step_at_time(5.0, 0.002, 0.5); }
    catch(const std::exception &error) { exception = error.what(); }

    const auto &d = controller.diagnostics();
    write_matrix(output / "H.csv", d.qpHessian);
    write_matrix(output / "f.csv", d.qpGradient);
    write_matrix(output / "B.csv", d.qpConstraintMatrix);
    write_matrix(output / "z.csv", d.qpConstraintVector);
    write_matrix(output / "U0.csv", d.qpSeed);
    write_matrix(output / "solution.csv", d.qpReturnedSolution);

    const Eigen::VectorXd slack = d.qpConstraintVector
        - d.qpConstraintMatrix * d.qpSeed;
    const double zeroThreshold = 1e-12;
    int seedBoundary = 0;
    int progressBoundary = 0;
    int contactBoundary = 0;
    const int baseVariables = 7 * static_cast<int>(horizon);
    const int baseRows = 2 * baseVariables + 1;
    for(int row = 0; row < slack.size(); ++row)
    {
        if(std::abs(slack(row)) > zeroThreshold) continue;
        ++seedBoundary;
        const bool progressUpper =
            row < baseVariables && row % 7 == 6;
        const bool progressLower =
            row >= baseVariables && row < 2 * baseVariables
            && (row - baseVariables) % 7 == 6;
        if(progressUpper || progressLower) ++progressBoundary;
        if(row >= baseRows) ++contactBoundary;
    }

    std::ofstream metadata(output / "metadata.json");
    metadata << std::setprecision(17)
             << "{\n"
             << "  \"horizon\": " << horizon << ",\n"
             << "  \"contact_mode\": \"" << modeName << "\",\n"
             << "  \"progress_configuration\": \""
             << (fixed ? "fixed_progress_h13" : "current_e07_variable_progress")
             << "\",\n"
             << "  \"progress_rate_min\": 0.0,\n"
             << "  \"progress_rate_nominal\": 0.1,\n"
             << "  \"progress_rate_max\": 0.125,\n"
             << "  \"seed_boundary_threshold\": " << zeroThreshold << ",\n"
             << "  \"seed_boundary_rows\": " << seedBoundary << ",\n"
             << "  \"seed_progress_boundary_rows\": " << progressBoundary << ",\n"
             << "  \"seed_contact_boundary_rows\": " << contactBoundary << ",\n"
             << "  \"converged\": " << (d.qpConverged ? "true" : "false") << ",\n"
             << "  \"iterations\": " << d.qpIterations << ",\n"
             << "  \"max_steps\": 50,\n"
             << "  \"final_step_size\": " << d.qpFinalStepSize << ",\n"
             << "  \"raw_primal_violation\": "
             << (d.qpConstraintMatrix * d.qpReturnedSolution
                 - d.qpConstraintVector).maxCoeff() << ",\n"
             << "  \"reported_primal_violation\": "
             << d.qpMaximumConstraintViolation << ",\n"
             << "  \"objective\": "
             << (0.5 * d.qpReturnedSolution.dot(
                     d.qpHessian * d.qpReturnedSolution)
                 + d.qpGradient.dot(d.qpReturnedSolution)) << ",\n"
             << "  \"runtime_exception\": \"" << exception << "\"\n"
             << "}\n";
    return d.qpHessian.size() == 0 ? 1 : 0;
}
