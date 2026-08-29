/**
 * @file SerialLinkLieAlgebraMPC.cpp
 * @brief Convex time-indexed MPC in SE(3) logarithmic error coordinates.
 */

#include <Control/TrajectoryTracking/SerialLinkLieAlgebraMPC.h>
#include "detail/LieAlgebraMPCPrediction.h"
#include <Math/CondensedMPC.h>
#include <Math/DiscreteIntegratorLQR.h>
#include <Math/MathFunctions.h>
#include <Math/QpAcceptance.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>

namespace RobotLibrary { namespace Control {

namespace {
Eigen::Matrix4d pose_matrix(RobotLibrary::Model::Pose pose)
{
    return pose.as_matrix();
}
}

SerialLinkLieAlgebraMPC::SerialLinkLieAlgebraMPC(
    std::shared_ptr<RobotLibrary::Model::KinematicTree> model,
    const std::string &endpointName,
    const RobotLibrary::Control::SerialLinkParameters &parameters,
    unsigned int horizon,
    double dt)
: SerialLinkTimeIndexedMPC(model, endpointName, parameters),
  _horizon(std::max(1u, horizon)),
  _dt(dt > 0.0 ? dt : 0.002),
  _qpSolver(parameters.qpsolver)
{
    _controlFrequency = 1.0 / _dt;
}

void
SerialLinkLieAlgebraMPC::set_feedback_bandwidth(const double positionGain,
                                                const double orientationGain)
{
    _wPosition = RobotLibrary::Math::integrator_stage_weight_for_gain(
        positionGain, _wLinearVelocity, _dt);
    _wOrientation = RobotLibrary::Math::integrator_stage_weight_for_gain(
        orientationGain, _wAngularVelocity, _dt);
    _terminalStateWeight.setZero();
    _terminalStateWeight.diagonal().head<3>().setConstant(
        RobotLibrary::Math::integrator_terminal_weight_for_gain(
            positionGain, _wLinearVelocity, _dt));
    _terminalStateWeight.diagonal().tail<3>().setConstant(
        RobotLibrary::Math::integrator_terminal_weight_for_gain(
            orientationGain, _wAngularVelocity, _dt));
    _useTerminalStateWeight = true;
    _warmStart.resize(0);
}

void
SerialLinkLieAlgebraMPC::set_trajectory(
    const RobotLibrary::Trajectory::CartesianSpline &trajectory)
{
    _trajectory = trajectory;
    _trajectorySet = true;
    _warmStart.resize(0);
}

void
SerialLinkLieAlgebraMPC::set_trajectory_frame(
    const RobotLibrary::Trajectory::CartesianTrajectoryFrameState &frame)
{
    RobotLibrary::Trajectory::validate_trajectory_frame(frame);
    _trajectoryFrame = frame;
}

void SerialLinkLieAlgebraMPC::clear_trajectory()
{
    _trajectorySet = false;
    _warmStart.resize(0);
}

Eigen::VectorXd
SerialLinkLieAlgebraMPC::track_endpoint_trajectory_at_time(const double &time)
{
    if(not _trajectorySet)
    {
        throw std::runtime_error(
            "[ERROR] [SERIAL LINK LIE ALGEBRA MPC] No trajectory has been set.");
    }

    update();
    const double start = _trajectory.start_time();
    const double end = _trajectory.end_time();
    const double t0 = std::clamp(std::isfinite(time) ? time : start, start, end);

    RobotLibrary::Trajectory::CartesianState currentReference =
        RobotLibrary::Trajectory::express_state_in_base(
            _trajectoryFrame, _trajectory.query_state(t0));
    RobotLibrary::Model::Pose currentPose = endpoint_pose();
    (void)pose_error(currentReference.pose);

    const Eigen::Matrix4d referenceTransform = pose_matrix(currentReference.pose);
    const Eigen::Matrix4d currentTransform = pose_matrix(currentPose);
    const LieAlgebraMPCVector initialError = RobotLibrary::Math::se3_logarithm(
        RobotLibrary::Math::se3_inverse(referenceTransform) * currentTransform);

    std::vector<LieAlgebraMPCVector> referenceBodyTwists;
    referenceBodyTwists.reserve(_horizon);
    for(unsigned int k = 0; k < _horizon; ++k)
    {
        const double sampleTime = std::clamp(t0 + static_cast<double>(k) * _dt, start, end);
        RobotLibrary::Trajectory::CartesianState reference =
            RobotLibrary::Trajectory::express_state_in_base(
                _trajectoryFrame, _trajectory.query_state(sampleTime));
        referenceBodyTwists.push_back(
            endpoint_twist_in_body(reference.pose.rotation(), reference.twist));
    }

    const auto solveStart = std::chrono::steady_clock::now();
    LieAlgebraMPCVector bodyTwist = solve_error_mpc(initialError, referenceBodyTwists);
    for(int i = 0; i < 3; ++i)
        bodyTwist(i) = std::clamp(bodyTwist(i), -_maxLinearSpeed, _maxLinearSpeed);
    for(int i = 3; i < 6; ++i)
        bodyTwist(i) = std::clamp(bodyTwist(i), -_maxAngularSpeed, _maxAngularSpeed);

    const LieAlgebraMPCVector baseTwist =
        endpoint_twist_in_base(currentPose.rotation(), bodyTwist);
    const Eigen::VectorXd jointCommand = resolve_endpoint_twist(baseTwist);
    record_time_indexed_diagnostics(
        baseTwist, bodyTwist, jointCommand, t0, _dt, _maxLinearSpeed, _maxAngularSpeed,
        std::chrono::duration<double>(std::chrono::steady_clock::now() - solveStart).count());
    return jointCommand;
}

Eigen::VectorXd
SerialLinkLieAlgebraMPC::track_endpoint_trajectory(
    const RobotLibrary::Model::Pose &desiredPose,
    const Eigen::Vector<double,6> &desiredVelocity,
    const Eigen::Vector<double,6> &desiredAcceleration)
{
    (void)desiredAcceleration;
    update();
    RobotLibrary::Model::Pose currentPose = endpoint_pose();
    RobotLibrary::Model::Pose desired = desiredPose;
    (void)pose_error(desiredPose);

    const LieAlgebraMPCVector initialError = RobotLibrary::Math::se3_logarithm(
        RobotLibrary::Math::se3_inverse(pose_matrix(desired)) * pose_matrix(currentPose));
    const LieAlgebraMPCVector referenceBodyTwist =
        endpoint_twist_in_body(desiredPose.rotation(), desiredVelocity);
    const std::vector<LieAlgebraMPCVector> referenceBodyTwists(
        _horizon, referenceBodyTwist);
    const LieAlgebraMPCVector bodyTwist =
        solve_error_mpc(initialError, referenceBodyTwists);
    return resolve_endpoint_twist(
        endpoint_twist_in_base(currentPose.rotation(), bodyTwist));
}

LieAlgebraMPCVector
SerialLinkLieAlgebraMPC::solve_error_mpc(
    const LieAlgebraMPCVector &initialError,
    const std::vector<LieAlgebraMPCVector> &referenceBodyTwists)
{
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    constexpr int dimension = 6;
    const int N = static_cast<int>(_horizon);
    const int stackedDimension = dimension * N;

    std::vector<MatrixXd> stateMatrices;
    std::vector<MatrixXd> inputMatrices;
    stateMatrices.reserve(_horizon);
    inputMatrices.reserve(_horizon);
    for(unsigned int k = 0; k < _horizon; ++k)
    {
        const LieAlgebraMPCVector reference =
            k < referenceBodyTwists.size()
            ? referenceBodyTwists[k]
            : LieAlgebraMPCVector::Zero();
        const auto [A, B] = lie_algebra_mpc_stage(reference, _dt);
        stateMatrices.emplace_back(A);
        inputMatrices.emplace_back(B);
    }
    const RobotLibrary::Math::CondensedPrediction prediction =
        RobotLibrary::Math::condense_time_varying_prediction(
            stateMatrices, inputMatrices);

    LieAlgebraMPCMatrix Q = LieAlgebraMPCMatrix::Zero();
    Q.diagonal().head<3>().setConstant(_wPosition);
    Q.diagonal().tail<3>().setConstant(_wOrientation);
    LieAlgebraMPCMatrix R = LieAlgebraMPCMatrix::Zero();
    R.diagonal().head<3>().setConstant(_wLinearVelocity);
    R.diagonal().tail<3>().setConstant(_wAngularVelocity);
    MatrixXd Qbar = RobotLibrary::Math::block_diagonal(Q, _horizon);
    const MatrixXd Rbar = RobotLibrary::Math::block_diagonal(R, _horizon);
    if(_useTerminalStateWeight)
    {
        Qbar.bottomRightCorner<dimension,dimension>() = _terminalStateWeight;
    }

    MatrixXd H = prediction.inputResponse.transpose() * Qbar
        * prediction.inputResponse + Rbar;
    H.diagonal().array() += 1e-8;
    const VectorXd f = prediction.inputResponse.transpose() * Qbar
        * prediction.stateTransition * initialError;

    VectorXd lower(stackedDimension);
    VectorXd upper(stackedDimension);
    for(int k = 0; k < N; ++k)
    {
        const LieAlgebraMPCVector reference =
            k < static_cast<int>(referenceBodyTwists.size())
            ? referenceBodyTwists[static_cast<size_t>(k)]
            : LieAlgebraMPCVector::Zero();
        lower.segment<3>(k * dimension).setConstant(-_maxLinearSpeed);
        upper.segment<3>(k * dimension).setConstant(_maxLinearSpeed);
        lower.segment<3>(k * dimension + 3).setConstant(-_maxAngularSpeed);
        upper.segment<3>(k * dimension + 3).setConstant(_maxAngularSpeed);
        lower.segment<dimension>(k * dimension) -= reference;
        upper.segment<dimension>(k * dimension) -= reference;
    }
    const RobotLibrary::Math::BoxConstraint box =
        RobotLibrary::Math::box_constraint(lower, upper);

    if(_warmStart.size() != stackedDimension)
        _warmStart = VectorXd::Zero(stackedDimension);
    for(int i = 0; i < stackedDimension; ++i)
        _warmStart(i) = std::clamp(_warmStart(i), lower(i), upper(i));
    VectorXd correction = _qpSolver.solve(
        H, f, box.constraintMatrix, box.constraintVector, _warmStart);
    const SolverResults<double> qpResults = _qpSolver.results();
    _timeIndexedDiagnostics.qpIterations = qpResults.numberOfSteps;
    _timeIndexedDiagnostics.qpFinalStepSize = qpResults.finalStepSize;
    _timeIndexedDiagnostics.qpConverged =
        qpResults.terminationReason == SolverTerminationReason::Converged;
    _timeIndexedDiagnostics.qpHitMaxIterations =
        qpResults.terminationReason == SolverTerminationReason::MaxIterations;
    const double violation =
        (box.constraintMatrix * correction - box.constraintVector).maxCoeff();
    _timeIndexedDiagnostics.qpPrimalViolation = std::max(0.0, violation);
    RobotLibrary::Math::require_qp_converged(
        qpResults, "[ERROR] [SERIAL LINK LIE ALGEBRA MPC] solve_error_mpc()");
    if(correction.size() != stackedDimension or not correction.allFinite())
    {
        throw std::runtime_error(
            "[ERROR] [SERIAL LINK LIE ALGEBRA MPC] QP returned an invalid solution.");
    }
    if(not std::isfinite(violation) or violation > 1e-6)
    {
        throw std::runtime_error(
            "[ERROR] [SERIAL LINK LIE ALGEBRA MPC] QP returned an infeasible solution.");
    }

    _warmStart = correction;
    const LieAlgebraMPCVector firstReference = referenceBodyTwists.empty()
        ? LieAlgebraMPCVector::Zero() : referenceBodyTwists.front();
    return firstReference + correction.head<dimension>();
}

} } // namespace RobotLibrary::Control
