/** E06-Q deterministic resolver and frozen-MPCC tolerance audit. */

#include <Control/Contact/SerialLinkPredictiveContactMPCC.h>
#include <Control/TrajectoryTracking/SerialLinkKinematic.h>
#include <Control/TrajectoryTracking/SerialLinkMPCC.h>
#include <Math/QPSolver.h>
#include <Model/KinematicTree.h>
#include <Trajectory/CartesianSpline.h>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using RobotLibrary::Control::MpccDiagnostics;
using RobotLibrary::Control::PredictiveContactMode;
using RobotLibrary::Control::PredictiveContactMpccDiagnostics;
using RobotLibrary::Control::PredictiveContactMpccParameters;
using RobotLibrary::Control::SerialLinkKinematic;
using RobotLibrary::Control::SerialLinkMPCC;
using RobotLibrary::Control::SerialLinkParameters;
using RobotLibrary::Control::SerialLinkPredictiveContactMPCC;
using RobotLibrary::Model::KinematicTree;
using RobotLibrary::Model::Pose;
using RobotLibrary::Trajectory::CartesianSpline;
using RobotLibrary::Trajectory::CartesianTrajectoryFrameState;

constexpr int kHorizon = 12;
constexpr double kDt = 0.002;
constexpr double kStiffness = 21303.75539503847;
const std::array<double, 8> kTolerances{
    5e-2, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8};
const std::array<double, 7> kLinearAmplitudes{
    1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3};

std::shared_ptr<KinematicTree> model_at(
    const std::string &urdf, const Eigen::VectorXd &q)
{
    auto model = std::make_shared<KinematicTree>(urdf);
    model->update_state(q, Eigen::VectorXd::Zero(q.size()));
    return model;
}

double manipulability(const Eigen::MatrixXd &jacobian)
{
    return std::sqrt(std::max(0.0, (jacobian * jacobian.transpose()).determinant()));
}

SerialLinkParameters parameters(const double tolerance)
{
    SerialLinkParameters result;
    result.controlFrequency = 500;
    result.minManipulability = 0.001;
    result.maxJointAcceleration = 10.0;
    result.qpsolver.method = "active set";
    result.qpsolver.stepSizeTolerance = tolerance;
    result.qpsolver.maxSteps = 100;
    return result;
}

Eigen::Matrix4d transform(const Eigen::Matrix3d &rotation,
                          const Eigen::Vector3d &translation)
{
    Eigen::Matrix4d result = Eigen::Matrix4d::Identity();
    result.block<3,3>(0,0) = rotation;
    result.block<3,1>(0,3) = translation;
    return result;
}

template<class Controller>
void set_path(Controller &controller, const Eigen::Matrix4d &parent)
{
    const Pose endpoint = controller.endpoint_pose();
    const Eigen::Matrix3d parentRotation = parent.block<3,3>(0,0);
    const Eigen::Vector3d parentPosition = parent.block<3,1>(0,3);
    const Pose start(
        parentRotation.transpose() * (endpoint.translation() - parentPosition),
        Eigen::Quaterniond(parentRotation.transpose()
                           * endpoint.quaternion().toRotationMatrix()));
    const Pose end(start.translation() + Eigen::Vector3d(0.03, 0.0, 0.0),
                   start.quaternion());
    controller.set_trajectory(CartesianSpline(
        start, end, Eigen::Vector<double,6>::Zero(), 0.0, 3.0));
    CartesianTrajectoryFrameState frame;
    frame.transformInBase = parent;
    controller.set_trajectory_frame(frame, 0.0);
}

void run_q0(const std::string &urdf, const std::filesystem::path &output)
{
    struct Configuration {const char *name; std::array<double,6> q;};
    const std::array<Configuration,3> configurations{{
        {"P0_contact", {1.5708, -1.5708, -1.5708, -3.14, -1.5708, 0.0}},
        {"P1_free", {0.0, -1.5708, 1.5708, -1.5708, -1.5708, 0.0}},
        {"P2_low_manip", {0.2, -1.30, 0.28, -1.72, -1.48, 0.15}}
    }};
    struct Direction {const char *name; Eigen::Vector<double,6> axis; bool rotation;};
    Eigen::Vector<double,6> normal = Eigen::Vector<double,6>::Zero();
    normal(1) = -1.0;
    Eigen::Vector<double,6> tangent = Eigen::Vector<double,6>::Zero();
    tangent(0) = 1.0;
    Eigen::Vector<double,6> rotation = Eigen::Vector<double,6>::Zero();
    rotation(3) = 1.0;
    const std::array<Direction,3> directions{{
        {"normal", normal, false}, {"tangent", tangent, false},
        {"rotation_x", rotation, true}}};

    std::ofstream stream(output);
    stream << std::setprecision(14)
           << "configuration,direction,requested_equivalent_m_s,requested_amplitude,"
              "tolerance,manipulability,singular_branch,resolved_directional,"
              "resolver_gain,relative_error,cross_leakage,qdot_norm,active_constraints,"
              "iterations,final_step_size,solve_time_us";
    for(int i = 0; i < 6; ++i) stream << ",qdot_" << i;
    stream << '\n';

    for(const auto &configuration : configurations) {
        const Eigen::VectorXd q = Eigen::Map<const Eigen::Vector<double,6>>(
            configuration.q.data());
        auto referenceModel = model_at(urdf, q);
        const auto frame = referenceModel->find_frame("wrist_3_link");
        const Eigen::MatrixXd jacobian = referenceModel->jacobian(frame);
        const double m = manipulability(jacobian);
        for(const auto &direction : directions) {
            for(const double linearAmplitude : kLinearAmplitudes) {
                const double requestedAmplitude = direction.rotation
                    ? linearAmplitude / 0.10 : linearAmplitude;
                const Eigen::Vector<double,6> request =
                    requestedAmplitude * direction.axis;
                for(const double tolerance : kTolerances) {
                    auto model = model_at(urdf, q);
                    SerialLinkKinematic controller(
                        model, "wrist_3_link", parameters(tolerance));
                    const auto start = std::chrono::steady_clock::now();
                    const Eigen::VectorXd qdot =
                        controller.resolve_endpoint_twist(request);
                    const double elapsed = 1e6 * std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - start).count();
                    const Eigen::Vector<double,6> realized = jacobian * qdot;
                    const double directional = direction.axis.dot(realized);
                    const double gain = directional / requestedAmplitude;
                    const double relativeError =
                        (realized - request).norm() / requestedAmplitude;
                    const double leakage =
                        (realized - directional * direction.axis).norm()
                        / requestedAmplitude;
                    Eigen::VectorXd lower(qdot.size()), upper(qdot.size());
                    // The requested amplitudes are far from physical bounds; count only
                    // numerically active box constraints from the URDF speed limits.
                    int active = 0;
                    for(int joint = 0; joint < qdot.size(); ++joint) {
                        const double speed = model->link(joint)->joint().speed_limit();
                        if(std::abs(std::abs(qdot(joint)) - speed) < 1e-9) ++active;
                    }
                    const auto solver = controller.results();
                    stream << configuration.name << ',' << direction.name << ','
                           << linearAmplitude << ',' << requestedAmplitude << ','
                           << tolerance << ',' << m << ',' << (m < 0.001) << ','
                           << directional << ',' << gain << ',' << relativeError << ','
                           << leakage << ',' << qdot.norm() << ',' << active << ','
                           << solver.numberOfSteps << ',' << solver.finalStepSize << ','
                           << elapsed;
                    for(int joint = 0; joint < qdot.size(); ++joint)
                        stream << ',' << qdot(joint);
                    stream << '\n';
                }
            }
        }
    }
}

struct Snapshot
{
    std::string name;
    MpccDiagnostics qp;
    PredictiveContactMpccDiagnostics contact;
    bool predictive{false};
};

Snapshot basic_snapshot(const std::string &urdf, const Eigen::VectorXd &q,
                        const bool translating)
{
    const Eigen::Matrix4d board = transform(
        Eigen::Matrix3d::Identity(), Eigen::Vector3d(-0.154, -0.454, 0.388));
    auto model = model_at(urdf, q);
    SerialLinkMPCC controller(
        model, "wrist_3_link", parameters(5e-2), kHorizon, kDt);
    set_path(controller, board);
    if(translating) {
        CartesianTrajectoryFrameState frame;
        frame.transformInBase = board;
        frame.transformInBase(1,3) -= 2.5e-5;
        controller.set_trajectory_frame(frame, kDt);
    }
    (void)controller.step_at_time(kDt, kDt);
    return {translating ? "S1_basic_translate" : "S0_basic_stationary",
            controller.diagnostics(), {}, false};
}

Snapshot predictive_snapshot(const std::string &urdf, const Eigen::VectorXd &q,
                             const bool pitch)
{
    const Eigen::Matrix4d board = transform(
        Eigen::Matrix3d::Identity(), Eigen::Vector3d(-0.154, -0.454, 0.388));
    auto model = model_at(urdf, q);
    SerialLinkPredictiveContactMPCC controller(
        model, "wrist_3_link", parameters(5e-2), kHorizon, kDt);
    PredictiveContactMpccParameters contactParameters;
    contactParameters.mode = PredictiveContactMode::Active;
    contactParameters.forceWeight = 1.0;
    controller.set_predictive_contact_parameters(contactParameters);
    controller.set_force_measurement(2.5, true);
    set_path(controller, board);
    if(pitch) {
        CartesianTrajectoryFrameState frame;
        frame.transformInBase = transform(
            Eigen::AngleAxisd(0.0015, Eigen::Vector3d::UnitX()).toRotationMatrix(),
            board.block<3,1>(0,3));
        controller.set_trajectory_frame(frame, kDt);
    }
    (void)controller.step_at_time(kDt, kDt);
    return {pitch ? "S3_predictive_pitch" : "S2_predictive_stationary",
            controller.diagnostics(), controller.predictive_contact_diagnostics(), true};
}

double percentile(std::vector<double> values, const double probability)
{
    std::sort(values.begin(), values.end());
    const size_t index = static_cast<size_t>(std::ceil(
        probability * static_cast<double>(values.size()))) - 1;
    return values[std::min(index, values.size() - 1)];
}

void run_q1(const std::string &urdf, const std::filesystem::path &output)
{
    const std::array<double,6> raw{
        1.5708, -1.5708, -1.5708, -3.14, -1.5708, 0.0};
    const Eigen::VectorXd q = Eigen::Map<const Eigen::Vector<double,6>>(raw.data());
    const std::array<Snapshot,3> snapshots{{
        basic_snapshot(urdf, q, false), basic_snapshot(urdf, q, true),
        predictive_snapshot(urdf, q, false)}};
    std::ofstream stream(output);
    stream << std::setprecision(14)
           << "snapshot,predictive,tolerance,dimension,constraint_rows,"
              "first_linear_error_m_s,first_angular_error_rad_s,progress_error,"
              "force_slack,penetration_slack_m,pred_force_1,pred_force_terminal,"
              "pred_force_min,pred_force_max,pred_delta_rho_terminal,objective,"
              "objective_gap,max_inequality_residual,iterations,final_step_size,"
              "solve_mean_us,solve_p99_us,solve_max_us";
    for(int i = 0; i < 7; ++i) stream << ",u0_" << i;
    stream << '\n';

    for(const Snapshot &snapshot : snapshots) {
        const auto &d = snapshot.qp;
        SolverOptions<double> referenceOptions;
        referenceOptions.method = "active set";
        referenceOptions.stepSizeTolerance = 1e-8;
        referenceOptions.maxSteps = 100;
        QPSolver<double> referenceSolver(referenceOptions);
        const Eigen::VectorXd reference = referenceSolver.solve(
            d.qpHessian, d.qpGradient, d.qpConstraintMatrix,
            d.qpConstraintVector, d.qpSeed);
        const double referenceObjective =
            reference.dot(0.5 * d.qpHessian * reference + d.qpGradient);
        for(const double tolerance : kTolerances) {
            SolverOptions<double> options;
            options.method = "active set";
            options.stepSizeTolerance = tolerance;
            options.maxSteps = 100;
            Eigen::VectorXd solution;
            SolverResults<double> solverResults;
            std::vector<double> timings;
            timings.reserve(300);
            for(int repetition = 0; repetition < 300; ++repetition) {
                QPSolver<double> solver(options);
                const auto start = std::chrono::steady_clock::now();
                solution = solver.solve(
                    d.qpHessian, d.qpGradient, d.qpConstraintMatrix,
                    d.qpConstraintVector, d.qpSeed);
                timings.push_back(1e6 * std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - start).count());
                solverResults = solver.results();
            }
            const Eigen::VectorXd difference = solution - reference;
            const double objective =
                solution.dot(0.5 * d.qpHessian * solution + d.qpGradient);
            const double violation = std::max(
                0.0, (d.qpConstraintMatrix * solution
                      - d.qpConstraintVector).maxCoeff());
            double forceSlack = 0.0;
            double penetrationSlack = 0.0;
            double force1 = 0.0, forceTerminal = 0.0;
            double forceMin = 0.0, forceMax = 0.0, rhoTerminal = 0.0;
            if(snapshot.predictive) {
                const int base = 7 * kHorizon;
                forceSlack = solution.segment(base, kHorizon).maxCoeff();
                penetrationSlack = 4e-6
                    * solution.segment(base + kHorizon, kHorizon).maxCoeff();
                const Eigen::VectorXd force = Eigen::VectorXd::Constant(kHorizon, 2.5)
                    + snapshot.contact.forceOffset
                    + snapshot.contact.forceMap * solution.head(base);
                force1 = force(0);
                forceTerminal = force(force.size() - 1);
                forceMin = force.minCoeff();
                forceMax = force.maxCoeff();
                rhoTerminal = (forceTerminal - 2.5) / kStiffness;
            }
            const double mean = std::accumulate(
                timings.begin(), timings.end(), 0.0) / timings.size();
            stream << snapshot.name << ',' << snapshot.predictive << ','
                   << tolerance << ',' << solution.size() << ','
                   << d.qpConstraintMatrix.rows() << ','
                   << difference.head<3>().norm() << ','
                   << difference.segment<3>(3).norm() << ','
                   << std::abs(difference(6)) << ',' << forceSlack << ','
                   << penetrationSlack << ',' << force1 << ',' << forceTerminal
                   << ',' << forceMin << ',' << forceMax << ',' << rhoTerminal
                   << ',' << objective << ',' << objective - referenceObjective
                   << ',' << violation << ',' << solverResults.numberOfSteps
                   << ',' << solverResults.finalStepSize << ',' << mean << ','
                   << percentile(timings, 0.99) << ','
                   << *std::max_element(timings.begin(), timings.end());
            for(int i = 0; i < 7; ++i) stream << ',' << solution(i);
            stream << '\n';
        }
    }
}

} // namespace

int main(int argc, char **argv)
{
    if(argc != 3) {
        std::cerr << "Usage: qp_resolution_audit URDF OUTPUT_DIRECTORY\n";
        return 2;
    }
    const std::filesystem::path output(argv[2]);
    std::filesystem::create_directories(output);
    run_q0(argv[1], output / "q0_resolver_sweep.csv");
    run_q1(argv[1], output / "q1_mpcc_sweep.csv");
    std::cout << "E06-Q deterministic audit written to " << output << '\n';
    return 0;
}
