/**
 * @file    SerialLinkMovingFrameMPC.cpp
 * @brief   Time-indexed MPC for trajectories expressed in a moving frame.
 */

#include <Control/SerialLinkMovingFrameMPC.h>
#include <Math/MathFunctions.h>
#include <Math/CondensedMPC.h>
#include <Math/QPSolver.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

using RobotLibrary::Math::quaternion_to_rotation_vector;

namespace RobotLibrary { namespace Control {

SerialLinkMovingFrameMPC::SerialLinkMovingFrameMPC(
    std::shared_ptr<RobotLibrary::Model::KinematicTree> model,
    const std::string &endpointName,
    const RobotLibrary::Control::SerialLinkParameters &parameters,
    unsigned int horizon,
    double dt)
: SerialLinkMPC(model, endpointName, parameters, horizon, dt)
{
}

RobotLibrary::Model::Pose
SerialLinkMovingFrameMPC::compose_reference_pose(
    const MovingFrameState &frame,
    const RobotLibrary::Model::Pose &relativePose)
{
    return frame.pose * relativePose;
}

Eigen::Vector<double,6>
SerialLinkMovingFrameMPC::compose_reference_twist(
    const MovingFrameState &frame,
    const RobotLibrary::Model::Pose &relativePose,
    const Eigen::Vector<double,6> &relativeTwist)
{
    const Eigen::Matrix3d frameRotation = frame.pose.rotation();
    const Eigen::Vector3d relativePositionBase = frameRotation * relativePose.translation();

    const Eigen::Vector3d frameLinearVelocity = frame.twist.head<3>();
    const Eigen::Vector3d frameAngularVelocity = frame.twist.tail<3>();
    const Eigen::Vector3d relativeLinearVelocityBase = frameRotation * relativeTwist.head<3>();
    const Eigen::Vector3d relativeAngularVelocityBase = frameRotation * relativeTwist.tail<3>();

    Eigen::Vector<double,6> twist = Eigen::Vector<double,6>::Zero();
    twist.head<3>() = frameLinearVelocity
                    + frameAngularVelocity.cross(relativePositionBase)
                    + relativeLinearVelocityBase;
    twist.tail<3>() = frameAngularVelocity + relativeAngularVelocityBase;
    return twist;
}

Eigen::VectorXd
SerialLinkMovingFrameMPC::track_moving_frame_trajectory_at_time(
    const double &time,
    const std::vector<MovingFrameState> &movingFramePrediction)
{
    if (!_trajectorySet)
    {
        throw std::runtime_error(
            "[ERROR] [SERIAL LINK MOVING FRAME MPC] track_moving_frame_trajectory_at_time(): "
            "No trajectory set. Call set_trajectory() first.");
    }

    if (movingFramePrediction.empty())
    {
        throw std::invalid_argument(
            "[ERROR] [SERIAL LINK MOVING FRAME MPC] track_moving_frame_trajectory_at_time(): "
            "Moving-frame prediction is empty.");
    }

    update();

    const unsigned int N = (_horizon == 0) ? 1 : _horizon;
    const double dt = (_dt > 0.0) ? _dt : 0.002;
    const double startTime = _trajectory.start_time();
    const double endTime = _trajectory.end_time();
    const double t0 = std::isfinite(time) ? std::clamp(time, startTime, endTime) : startTime;
    const auto predictionAt = [&](const unsigned int index) -> const MovingFrameState&
    {
        return movingFramePrediction[std::min<std::size_t>(index, movingFramePrediction.size() - 1)];
    };

    RobotLibrary::Model::Pose currentPose = endpoint_pose();

    Eigen::Matrix<double,6,1> x0;
    x0.head<3>() = currentPose.translation();
    x0.tail<3>() = quaternion_to_rotation_vector(currentPose.quaternion());

    std::vector<Eigen::Matrix<double,6,1>> xRefStack;
    std::vector<Eigen::Matrix<double,6,1>> uRefStack;
    xRefStack.reserve(N);
    uRefStack.reserve(N);

    for (unsigned int k = 0; k < N; ++k)
    {
        const double stateTime = std::clamp(t0 + static_cast<double>(k + 1) * dt, startTime, endTime);
        const double controlTime = std::clamp(t0 + static_cast<double>(k) * dt, startTime, endTime);

        RobotLibrary::Trajectory::CartesianState stateReference = _trajectory.query_state(stateTime);
        RobotLibrary::Trajectory::CartesianState controlReference = _trajectory.query_state(controlTime);

        const MovingFrameState &stateFrame = predictionAt(k + 1);
        const MovingFrameState &controlFrame = predictionAt(k);
        const RobotLibrary::Model::Pose referencePoseBase =
            compose_reference_pose(stateFrame, stateReference.pose);

        Eigen::Matrix<double,6,1> xRef;
        xRef.head<3>() = referencePoseBase.translation();
        xRef.tail<3>() = quaternion_to_rotation_vector(referencePoseBase.quaternion());

        xRefStack.push_back(xRef);
        uRefStack.push_back(compose_reference_twist(controlFrame,
                                                    controlReference.pose,
                                                    controlReference.twist));
    }

    Eigen::Matrix<double,6,1> uOpt =
        (_contactParameters.mode == ContactMode::Disabled)
        ? solveMPC(x0, xRefStack, uRefStack)
        : solve_contact_mpc(x0, xRefStack, uRefStack, movingFramePrediction);

    for (int i = 0; i < 3; ++i)
    {
        uOpt[i] = std::clamp(uOpt[i], -_maxLinearSpeed, _maxLinearSpeed);
    }
    for (int i = 3; i < 6; ++i)
    {
        uOpt[i] = std::clamp(uOpt[i], -_maxAngularSpeed, _maxAngularSpeed);
    }

    return resolve_endpoint_twist(uOpt);
}

  ///////////////////////////////////////////////////////////////////////////////////////////////////
 //                          Track using the in-class board predictor                             //
///////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::VectorXd
SerialLinkMovingFrameMPC::track_moving_frame_trajectory_at_time(const double &time)
{
    if (_boardSamples.empty())
    {
        throw std::runtime_error(
            "[ERROR] [SERIAL LINK MOVING FRAME MPC] track_moving_frame_trajectory_at_time(): "
            "No board pose samples. Call update_board_pose() first.");
    }

    const unsigned int N = (_horizon == 0) ? 1 : _horizon;

    return track_moving_frame_trajectory_at_time(time, predict_board_motion(N));
}

  ///////////////////////////////////////////////////////////////////////////////////////////////////
 //                       Ingest a board pose sample; estimate v, a, omega                        //
///////////////////////////////////////////////////////////////////////////////////////////////////
void
SerialLinkMovingFrameMPC::update_board_pose(const RobotLibrary::Model::Pose &boardPose, double time)
{
    _boardSamples.push_back({boardPose, time});
    while (_boardSamples.size() > 3) _boardSamples.pop_front();

    const std::size_t m = _boardSamples.size();

    if (m < 2)
    {
        _boardLinearVelocity.setZero();
        _boardLinearAcceleration.setZero();
        _boardAngularVelocity.setZero();
        _boardEstimateValid = (m >= 1);
        return;
    }

    const BoardSample &last = _boardSamples[m - 1];
    const BoardSample &prev = _boardSamples[m - 2];

    const double dt1 = last.time - prev.time;
    if (dt1 < 1e-6)
    {
        _boardEstimateValid = true;                                                                 // Stale/duplicate stamp: keep previous estimate
        return;
    }

    _boardLinearVelocity = (last.pose.translation() - prev.pose.translation()) / dt1;

    const Eigen::Quaterniond qRelative = last.pose.quaternion() * prev.pose.quaternion().inverse();
    _boardAngularVelocity = RobotLibrary::Math::quaternion_to_rotation_vector(qRelative) / dt1;

    _boardLinearAcceleration.setZero();

    if (m >= 3 and _contactParameters.useAcceleration)
    {
        const BoardSample &prev2 = _boardSamples[m - 3];
        const double dt0 = prev.time - prev2.time;

        if (dt0 > 1e-6)
        {
            const Eigen::Vector3d olderVelocity =
                (prev.pose.translation() - prev2.pose.translation()) / dt0;

            Eigen::Vector3d acceleration =
                (_boardLinearVelocity - olderVelocity) / (0.5 * (dt0 + dt1));

            const double maxA = _contactParameters.maxBoardAcceleration;
            for (int i = 0; i < 3; ++i) acceleration(i) = std::clamp(acceleration(i), -maxA, maxA);

            _boardLinearAcceleration = acceleration;
        }
    }

    _boardEstimateValid = true;
}

  ///////////////////////////////////////////////////////////////////////////////////////////////////
 //                  Roll out the board motion (constant acceleration model)                      //
///////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<MovingFrameState>
SerialLinkMovingFrameMPC::predict_board_motion(unsigned int steps) const
{
    std::vector<MovingFrameState> prediction;
    prediction.reserve(static_cast<std::size_t>(steps) + 1);

    if (_boardSamples.empty()) return prediction;

    const double dt = (_dt > 0.0) ? _dt : 0.002;

    const RobotLibrary::Model::Pose &latest = _boardSamples.back().pose;
    const Eigen::Vector3d    p0 = latest.translation();
    const Eigen::Quaterniond q0 = latest.quaternion();

    const Eigen::Vector3d v = _boardLinearVelocity;
    const Eigen::Vector3d a = _boardLinearAcceleration;
    const Eigen::Vector3d w = _boardAngularVelocity;

    for (unsigned int k = 0; k <= steps; ++k)
    {
        const double tau = static_cast<double>(k) * dt;

        MovingFrameState state;

        const Eigen::Matrix3d Rk = RobotLibrary::Math::so3_exponential(w * tau) * q0.toRotationMatrix();

        state.pose = RobotLibrary::Model::Pose(p0 + v * tau + 0.5 * a * tau * tau,
                                               Eigen::Quaterniond(Rk).normalized());

        state.twist.head<3>() = v + a * tau;
        state.twist.tail<3>() = w;

        prediction.push_back(state);
    }

    return prediction;
}

  ///////////////////////////////////////////////////////////////////////////////////////////////////
 //         Force-aware MPC: surface-frame tracking + board-relative normal-force term            //
///////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Matrix<double,6,1>
SerialLinkMovingFrameMPC::solve_contact_mpc(const Eigen::Matrix<double,6,1> &x0,
                                            const std::vector<Eigen::Matrix<double,6,1>> &xRefStack,
                                            const std::vector<Eigen::Matrix<double,6,1>> &uRefStack,
                                            const std::vector<MovingFrameState> &boardStack)
{
    using namespace Eigen;

    const ContactParameters &cp = _contactParameters;
    const bool useConstraint = (cp.mode == ContactMode::Constraint);

    const int nu = 6;
    const unsigned int N  = (_horizon == 0) ? 1 : _horizon;
    const double       dt = (_dt > 0.0) ? _dt : 0.002;

    const int twistVars      = nu * static_cast<int>(N);
    const int forceSlackVars = useConstraint ? static_cast<int>(N) : 0;
    const int V              = twistVars + forceSlackVars;

    const Eigen::Vector<double,6> uRefFallback =
        uRefStack.empty() ? Eigen::Vector<double,6>::Zero() : uRefStack.back();
    auto uRefAt = [&](unsigned int k) -> Eigen::Vector<double,6>
    {
        return (k < uRefStack.size()) ? uRefStack[k] : uRefFallback;
    };
    auto boardAt = [&](unsigned int idx) -> const MovingFrameState&
    {
        return boardStack[std::min<std::size_t>(idx, boardStack.size() - 1)];
    };

    // Condensed tool prediction X = Ax * x0 + Bu * U with A = I, B = dt I (state = [pos; rotvec]).
    MatrixXd A = MatrixXd::Identity(6, 6);
    MatrixXd B = dt * MatrixXd::Identity(6, nu);
    const MatrixXd &Bu = RobotLibrary::Math::condense_prediction(A, B, N).inputResponse;

    // Surface frame from the board orientation at the current instant.
    const Matrix3d Rboard0 = boardAt(0).pose.rotation();
    const Vector3d outwardNormal = (Rboard0 * cp.normalAxisInBoard).normalized();
    const Vector3d inwardNormal  = -outwardNormal;

    Vector3d tangentX = Rboard0 * cp.tangentAxisInBoard;
    tangentX -= inwardNormal * inwardNormal.dot(tangentX);
    if (tangentX.norm() < 1e-6)
    {
        tangentX = Vector3d::UnitX() - inwardNormal * inwardNormal.dot(Vector3d::UnitX());
        if (tangentX.norm() < 1e-6) tangentX = Vector3d::UnitY() - inwardNormal * inwardNormal.dot(Vector3d::UnitY());
    }
    tangentX.normalize();
    const Vector3d tangentY = inwardNormal.cross(tangentX).normalized();

    Matrix3d surfaceBasis;
    surfaceBasis.row(0) = tangentX.transpose();
    surfaceBasis.row(1) = tangentY.transpose();
    surfaceBasis.row(2) = inwardNormal.transpose();

    const Vector3d toolP0  = x0.head<3>();
    const Vector3d boardP0 = boardAt(0).pose.translation();

    MatrixXd H = MatrixXd::Zero(V, V);
    VectorXd g = VectorXd::Zero(V);

    // ---- Pose tracking (surface-weighted) + force loss --------------------------------------------
    for (unsigned int k = 0; k < N; ++k)
    {
        MatrixXd posMap = MatrixXd::Zero(3, V);                                                      // tool position vs U
        posMap.block(0, 0, 3, twistVars) = Bu.block(6 * k, 0, 3, twistVars);
        MatrixXd oriMap = MatrixXd::Zero(3, V);                                                      // tool rotvec vs U
        oriMap.block(0, 0, 3, twistVars) = Bu.block(6 * k + 3, 0, 3, twistVars);

        // Position error e = refPos - toolPos = (refPos - toolP0) - posMap*U, rotated to surface frame.
        const MatrixXd surfPosMap   = surfaceBasis * (-posMap);
        const Vector3d surfPosConst = surfaceBasis * (xRefStack[k].head<3>() - toolP0);
        Matrix3d Wpos = Matrix3d::Zero();
        Wpos(0, 0) = cp.tangentPositionWeight;
        Wpos(1, 1) = cp.tangentPositionWeight;
        Wpos(2, 2) = cp.normalPositionWeight;
        H += 2.0 * surfPosMap.transpose() * Wpos * surfPosMap;
        g += 2.0 * surfPosMap.transpose() * Wpos * surfPosConst;

        // Orientation error (absolute rotation-vector difference, as in the parent solveMPC).
        const MatrixXd oriErrMap   = -oriMap;
        const Vector3d oriErrConst = xRefStack[k].tail<3>() - x0.tail<3>();
        const Matrix3d Wori = cp.orientationWeight * Matrix3d::Identity();
        H += 2.0 * oriErrMap.transpose() * Wori * oriErrMap;
        g += 2.0 * oriErrMap.transpose() * Wori * oriErrConst;

        // Board-relative normal force: F_k = F_meas + k_c*( n^T*toolDisp(U) - n^T*boardDisp_k ).
        if (cp.mode == ContactMode::Loss)
        {
            const RowVectorXd forceMap = cp.forceResponseGain * inwardNormal.transpose() * posMap;
            const double boardDisp = inwardNormal.dot(boardAt(k + 1).pose.translation() - boardP0);
            const double forceConst = _measuredNormalForce - cp.forceResponseGain * boardDisp;
            const double residualConst = forceConst - cp.targetForce;
            H += 2.0 * cp.forceWeight * forceMap.transpose() * forceMap;
            g += 2.0 * cp.forceWeight * forceMap.transpose() * residualConst;
        }
        else if (useConstraint)
        {
            const int slackIdx = twistVars + static_cast<int>(k);
            H(slackIdx, slackIdx) += 2.0 * cp.slackWeight;
        }
    }

    // ---- Velocity + delta-velocity regularisation -------------------------------------------------
    const Eigen::Vector<double,6> prevCommand =
        _hasPreviousContactCommand ? _previousContactCommand : uRefAt(0);

    for (unsigned int k = 0; k < N; ++k)
    {
        const int idx = 6 * static_cast<int>(k);
        const Matrix<double,6,6> Wv = cp.velocityWeight * Matrix<double,6,6>::Identity();
        H.block(idx, idx, 6, 6) += 2.0 * Wv;
        g.segment(idx, 6)      += -2.0 * Wv * uRefAt(k);

        MatrixXd deltaMap = MatrixXd::Zero(6, V);
        Eigen::Vector<double,6> deltaConst = Eigen::Vector<double,6>::Zero();
        deltaMap.block(0, idx, 6, 6) = Matrix<double,6,6>::Identity();
        if (k == 0) deltaConst = -prevCommand;
        else        deltaMap.block(0, idx - 6, 6, 6) = -Matrix<double,6,6>::Identity();
        H += 2.0 * cp.deltaVelocityWeight * deltaMap.transpose() * deltaMap;
        g += 2.0 * cp.deltaVelocityWeight * deltaMap.transpose() * deltaConst;
    }

    H += 1e-8 * MatrixXd::Identity(V, V);

    // ---- Inequalities: velocity box (+ force band with slack) -------------------------------------
    std::vector<RowVectorXd> rows;
    std::vector<double>      bounds;
    auto addInequality = [&](const RowVectorXd &row, double bound) { rows.push_back(row); bounds.push_back(bound); };

    for (unsigned int k = 0; k < N; ++k)
    {
        for (int axis = 0; axis < 6; ++axis)
        {
            const double limit = (axis < 3) ? _maxLinearSpeed : _maxAngularSpeed;
            RowVectorXd e = RowVectorXd::Zero(V);
            e(6 * static_cast<int>(k) + axis) = 1.0;
            addInequality(e, limit);
            addInequality(-e, limit);
        }

        if (useConstraint)
        {
            MatrixXd posMap = MatrixXd::Zero(3, V);
            posMap.block(0, 0, 3, twistVars) = Bu.block(6 * k, 0, 3, twistVars);
            const RowVectorXd forceMap = cp.forceResponseGain * inwardNormal.transpose() * posMap;
            const double boardDisp  = inwardNormal.dot(boardAt(k + 1).pose.translation() - boardP0);
            const double forceConst = _measuredNormalForce - cp.forceResponseGain * boardDisp;
            const int    slackIdx   = twistVars + static_cast<int>(k);

            RowVectorXd upper = forceMap; upper(slackIdx) = -1.0;                                    // F_k - s <= F_d + tol
            addInequality(upper, cp.targetForce + cp.forceTolerance - forceConst);

            RowVectorXd lower = -forceMap; lower(slackIdx) = -1.0;                                   // -(F_k) - s <= -(F_d - tol)
            addInequality(lower, forceConst - (cp.targetForce - cp.forceTolerance));

            RowVectorXd slackLow = RowVectorXd::Zero(V); slackLow(slackIdx) = -1.0;                  // s >= 0
            addInequality(slackLow, 0.0);

            RowVectorXd slackHigh = RowVectorXd::Zero(V); slackHigh(slackIdx) = 1.0;                 // s <= maxSlack
            addInequality(slackHigh, cp.maxForceSlack);
        }
    }

    MatrixXd inequalities(static_cast<int>(rows.size()), V);
    VectorXd bound(static_cast<int>(bounds.size()));
    for (std::size_t i = 0; i < rows.size(); ++i)
    {
        inequalities.row(static_cast<int>(i)) = rows[i];
        bound(static_cast<int>(i))            = bounds[i];
    }

    // Warm start.
    if (_contactWarmStart.size() != V)
    {
        _contactWarmStart = VectorXd::Zero(V);
        for (unsigned int k = 0; k < N; ++k) _contactWarmStart.segment(6 * k, 6) = uRefAt(k);
    }

    Eigen::Vector<double,6> command = uRefAt(0);

    try
    {
        SolverOptions<double> options;
        options.method   = "active set";
        options.maxSteps = 10;
        QPSolver<double> solver(options);

        const VectorXd solution = solver.solve(H, g, inequalities, bound, _contactWarmStart);

        if (solution.size() == V)
        {
            _contactWarmStart = solution;
            command = solution.head<6>();
        }
    }
    catch (const std::exception &)
    {
        // Fall back to the board-glued feedforward twist on solver failure.
    }

    _previousContactCommand    = command;
    _hasPreviousContactCommand = true;

    return command;
}

} } // namespace
