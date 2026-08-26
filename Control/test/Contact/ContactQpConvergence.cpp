/** @file ContactQpConvergence.cpp
 *
 * Deployment-horizon convergence contract for the contact MPCC QP.
 *
 * Every other Control/test/Contact case builds its controller at horizon 1 to 4,
 * where the condensed QP is small enough to converge inside a handful of
 * active-set iterations. A solver-conditioning defect is a function of problem
 * size, so those cases cannot see one. This file builds at the deployment
 * horizon and asserts on the solver's termination state rather than only on the
 * direction of the resulting command: a truncated solve still returns a
 * feasible vector that passes every direction and finiteness check.
 */

#include <Control/Contact/SerialLinkContactMPCC.h>
#include <Math/QPSolver.h>
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

constexpr int kHorizon = 20;
constexpr double kDt = 0.002;
constexpr unsigned int kDeploymentMaxSteps = 50;   // ROS control_parameters.yaml
constexpr double kStepSizeTolerance = 1e-5;        // ROS control_parameters.yaml

int failures = 0;

void check(const bool condition, const std::string &message)
{
    if(not condition)
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
            / ("contact_qp_convergence_" + std::to_string(nonce) + ".urdf");
        std::ofstream stream(path);
        stream << R"(<robot name="contact_qp_convergence_fixture">
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

SerialLinkParameters deployment_parameters()
{
    SerialLinkParameters parameters;
    parameters.qpsolver.maxSteps = kDeploymentMaxSteps;
    parameters.qpsolver.stepSizeTolerance = kStepSizeTolerance;
    return parameters;
}

ContactParameters contact_parameters(const ContactMode mode)
{
    ContactParameters parameters;
    parameters.mode = mode;
    parameters.normalAxisInBoard = Eigen::Vector3d::UnitZ();
    parameters.tangentAxisInBoard = Eigen::Vector3d::UnitX();
    parameters.targetForce = 5.0;
    parameters.forceResponseGain = 1500.0;
    return parameters;
}

void configure(SerialLinkContactMPCC &controller, const ContactMode mode)
{
    const Pose start = controller.endpoint_pose();
    const Pose finish(start.translation() + Eigen::Vector3d(0.15, 0.0, 0.0),
                      start.quaternion());
    controller.set_trajectory(
        CartesianSpline(start, finish, Eigen::Vector<double,6>::Zero(), 0.0, 1.0));
    controller.set_contact_parameters(contact_parameters(mode));
    ContactState state;
    state.inContact = true;
    state.normalBase = Eigen::Vector3d::UnitZ();
    state.measuredNormalForce = 3.0;
    controller.set_contact_state(state);
}

/**
 * The deployment-horizon QP must reach its convergence exit inside the
 * configured iteration budget, and must say so.
 */
void converges_at_deployment_horizon(const UrdfFixture &fixture,
                                     const ContactMode mode,
                                     const std::string &label)
{
    auto model = make_model(fixture.path);
    SerialLinkContactMPCC controller(
        model, "tool", deployment_parameters(), kHorizon, kDt);
    configure(controller, mode);

    for(int cycle = 0; cycle < 25; ++cycle)
    {
        (void)controller.step(kDt);
        const auto &diagnostics = controller.diagnostics();
        check(diagnostics.qpConverged,
              label + ": QP converged within the deployment iteration budget");
        check(diagnostics.qpIterations < static_cast<double>(kDeploymentMaxSteps),
              label + ": QP did not exhaust its iteration budget");
        check(diagnostics.qpFinalStepSize <= kStepSizeTolerance,
              label + ": final step size is inside the configured tolerance");
        check(controller.contact_diagnostics().qpConverged,
              label + ": contact diagnostics report the same convergence state");
        if(failures > 0) return;                                                                    // Keep the first failure legible
    }
}

/**
 * Every stage's pinned progress rate must be carried by exactly one equality
 * row, not by two opposing inequalities.
 */
void pinned_progress_becomes_equalities(const UrdfFixture &fixture)
{
    auto model = make_model(fixture.path);
    SerialLinkContactMPCC controller(
        model, "tool", deployment_parameters(), kHorizon, kDt);
    configure(controller, ContactMode::Loss);
    (void)controller.step(kDt);

    const auto &diagnostics = controller.diagnostics();
    check(diagnostics.qpEqualityRows == static_cast<double>(kHorizon),
          "one equality row per pinned progress rate");
    check(diagnostics.qpEqualityViolation <= 1e-9,
          "the returned solution honours every fixed-variable equality");

    // No retained inequality row may still pin a variable between coincident
    // bounds; that is the degenerate pair the split exists to remove.
    const Eigen::MatrixXd &B = diagnostics.qpConstraintMatrix;
    int coincidentPairs = 0;
    for(Eigen::Index i = 0; i < B.rows(); ++i)
    {
        for(Eigen::Index j = i + 1; j < B.rows(); ++j)
        {
            if((B.row(i) + B.row(j)).cwiseAbs().maxCoeff() > 1e-14) continue;
            if(std::abs(diagnostics.qpConstraintVector(i)
                        + diagnostics.qpConstraintVector(j)) <= 1e-14)
            {
                ++coincidentPairs;
            }
        }
    }
    check(coincidentPairs == 0,
          "no opposing inequality pair still fixes a variable");
}

/**
 * A variable that still has a genuine interval must keep both bounds. Near the
 * end of the path the progress lower bound relaxes below the upper one, and
 * pinning it there would replace a conditioning defect with a wrong answer.
 */
void relaxed_progress_stays_free(const UrdfFixture &fixture)
{
    auto model = make_model(fixture.path);
    SerialLinkContactMPCC controller(
        model, "tool", deployment_parameters(), kHorizon, kDt);
    configure(controller, ContactMode::Loss);
    (void)controller.step(kDt, 0.999);

    check(controller.diagnostics().qpEqualityRows == 0.0,
          "a genuinely bounded progress rate is not pinned near the path end");
}

/**
 * The split must not change the answer. Restoring each equality as the opposing
 * inequality pair it replaced, and solving that original formulation to
 * convergence, must reproduce the same optimum.
 */
void split_preserves_the_optimum(const UrdfFixture &fixture, const ContactMode mode,
                                 const std::string &label)
{
    auto model = make_model(fixture.path);
    SerialLinkContactMPCC controller(
        model, "tool", deployment_parameters(), kHorizon, kDt);
    configure(controller, mode);
    (void)controller.step(kDt);

    const auto &d = controller.diagnostics();
    const Eigen::Index equalityRows = d.qpEqualityMatrix.rows();
    check(equalityRows > 0, label + ": the frozen QP carries an equality block");
    if(equalityRows == 0) return;

    Eigen::MatrixXd B(d.qpConstraintMatrix.rows() + 2 * equalityRows,
                      d.qpConstraintMatrix.cols());
    Eigen::VectorXd z(B.rows());
    B.topRows(d.qpConstraintMatrix.rows()) = d.qpConstraintMatrix;
    z.head(d.qpConstraintMatrix.rows()) = d.qpConstraintVector;
    B.middleRows(d.qpConstraintMatrix.rows(), equalityRows) = d.qpEqualityMatrix;
    z.segment(d.qpConstraintMatrix.rows(), equalityRows) = d.qpEqualityVector;
    B.bottomRows(equalityRows) = -d.qpEqualityMatrix;
    z.tail(equalityRows) = -d.qpEqualityVector;

    SolverOptions<double> options;
    options.method = "active set";
    options.stepSizeTolerance = kStepSizeTolerance;
    options.maxSteps = 4000;                                                                        // Whatever the original formulation needs
    QPSolver<double> reference(options);
    const Eigen::VectorXd original =
        reference.solve(d.qpHessian, d.qpGradient, B, z, d.qpSeed);

    check(reference.results().converged,
          label + ": the restored original formulation converges when given enough steps");
    check((original - d.qpReturnedSolution).cwiseAbs().maxCoeff() <= 1e-9,
          label + ": the split reaches the same optimum as the original box");
}

} // namespace

int main()
{
    const UrdfFixture fixture;

    converges_at_deployment_horizon(fixture, ContactMode::Loss, "Loss");
    converges_at_deployment_horizon(fixture, ContactMode::Constraint, "Constraint");
    pinned_progress_becomes_equalities(fixture);
    relaxed_progress_stays_free(fixture);
    split_preserves_the_optimum(fixture, ContactMode::Loss, "Loss");
    split_preserves_the_optimum(fixture, ContactMode::Constraint, "Constraint");

    if(failures == 0) std::cout << "contact_qp_convergence_test: all checks passed\n";
    return failures == 0 ? 0 : 1;
}
