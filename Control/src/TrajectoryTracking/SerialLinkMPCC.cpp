/**
 * @file    SerialLinkMPCC.cpp
 * @brief   MPCC controller implementation for serial link robot arms.
 */

#include <Control/TrajectoryTracking/SerialLinkMPCC.h>
#include "detail/ProgressSchedule.h"
#include "detail/RmpccCostGeometry.h"
#include "detail/RealtimeRecoveryQpCapture.h"
#include <Math/CondensedMPC.h>
#include <Math/BoxAwareActiveSet.h>
#include <Math/DiscreteIntegratorLQR.h>
#include <Math/MathFunctions.h>

#include <Eigen/Geometry>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <vector>

using RobotLibrary::Math::quaternion_to_rotation_vector;

namespace RobotLibrary { namespace Control {

SerialLinkMPCC::SerialLinkMPCC(std::shared_ptr<RobotLibrary::Model::KinematicTree> model,
                               const std::string &endpointName,
                               const RobotLibrary::Control::SerialLinkParameters &parameters,
                               unsigned int horizon,
                               double dt)
: SerialLinkVelocityBase(model, endpointName, parameters),
  _horizon((horizon == 0) ? 1 : horizon),
  _dt((dt > 0.0) ? dt : 1.0 / std::max(parameters.controlFrequency, 1u)),
  _qpSolver(parameters.qpsolver),
  _qpOptions(parameters.qpsolver)
{
    // The action-loop frequency is shared by all velocity controllers. The
    // prediction step follows it by default, but an explicit dt may still be
    // supplied for library-only use without silently changing the loop rate.
    _matchedPositionGain = parameters.cartesianPoseGain.diagonal().head<3>().mean();
    _matchedOrientationGain =
        0.5 * parameters.cartesianPoseGain.diagonal().tail<3>().mean();
}

Eigen::VectorXd
SerialLinkMPCC::track_endpoint_trajectory(const RobotLibrary::Model::Pose &desiredPose,
                                          const Eigen::Vector<double,6> &desiredVelocity,
                                          const Eigen::Vector<double,6> &desiredAcceleration)
{
    (void)desiredPose;
    (void)desiredVelocity;
    (void)desiredAcceleration;

    throw std::logic_error(
        "[ERROR] [SERIAL LINK MPCC] track_endpoint_trajectory(): "
        "Single-pose tracking is not MPCC. Call set_trajectory() and then step(dt).");
}

void
SerialLinkMPCC::set_trajectory(const RobotLibrary::Trajectory::CartesianSpline &trajectory)
{
    const double duration = trajectory.end_time() - trajectory.start_time();
    if(not std::isfinite(duration) or duration <= 0.0)
    {
        throw std::invalid_argument(
            "[ERROR] [SERIAL LINK MPCC] set_trajectory(): Trajectory duration must be positive.");
    }

    _trajectory = trajectory;
    _trajectorySet = true;
    _pathProgress = 0.0;
    _uLast.setZero();
    _warmStart.resize(0);
    _parentFrameMotion.reset();
    _diagnostics = MpccDiagnostics();
    _vProgressNominal = 1.0 / duration;
    _vProgressMin = (_ablationProfile == MpccAblationProfile::Baseline
                     or _ablationProfile == MpccAblationProfile::MatchedFeedbackGain)
        ? _vProgressNominal
        : (_ablationProfile == MpccAblationProfile::OptimizedProgress)
        ? 0.25 * _vProgressNominal : 0.0;
    if(_ablationProfile == MpccAblationProfile::FeedforwardCorrection
       or _ablationProfile == MpccAblationProfile::MatchedFeedbackGain)
    {
        _vProgressMax = _vProgressNominal;
    }
    else if(_ablationProfile == MpccAblationProfile::OptimizedProgress)
    {
        // The supplied trajectory remains the fastest admissible schedule.
        // MPCC may slow down for contouring accuracy, but may not silently
        // turn a time-indexed reference into a faster trajectory.
        _vProgressMax = _vProgressNominal;
    }
    else
    {
        _vProgressMax = std::max(0.1, _vProgressNominal);
    }
}

void
SerialLinkMPCC::set_ablation_profile(const MpccAblationProfile profile)
{
    if(_trajectorySet)
    {
        throw std::logic_error(
            "[ERROR] [SERIAL LINK MPCC] set_ablation_profile(): "
            "Configure the profile before set_trajectory().");
    }
    _ablationProfile = profile;
}

void
SerialLinkMPCC::set_angular_velocity_limit(const double angularVelocityMax)
{
    if(not std::isfinite(angularVelocityMax) or angularVelocityMax <= 0.0)
    {
        throw std::invalid_argument(
            "[ERROR] [SERIAL LINK MPCC] set_angular_velocity_limit(): "
            "Limit must be finite and positive.");
    }
    _vMaxAngular = angularVelocityMax;
}

void
SerialLinkMPCC::set_linear_velocity_limit(const double linearVelocityMax)
{
    if(not std::isfinite(linearVelocityMax) or linearVelocityMax <= 0.0)
    {
        throw std::invalid_argument(
            "[ERROR] [SERIAL LINK MPCC] set_linear_velocity_limit(): "
            "Limit must be finite and positive.");
    }
    _vMaxLinear = linearVelocityMax;
}

void
SerialLinkMPCC::set_trajectory_frame(
    const RobotLibrary::Trajectory::CartesianTrajectoryFrameState &frame)
{
    RobotLibrary::Trajectory::validate_trajectory_frame(frame);
    _trajectoryFrame = frame;
    _parentFrameMotion.set_static_pose(frame.transformInBase);
}

void
SerialLinkMPCC::set_trajectory_frame(
    const RobotLibrary::Trajectory::CartesianTrajectoryFrameState &frame,
    const double timestampSeconds)
{
    RobotLibrary::Trajectory::validate_trajectory_frame(frame);
    _trajectoryFrame = frame;
    _parentFrameMotion.update(frame.transformInBase, timestampSeconds);
}

Eigen::VectorXd
SerialLinkMPCC::step(const double dt)
{
    using Clock = std::chrono::steady_clock;
    const auto totalStart = _timingDiagnosticsEnabled ? Clock::now() : Clock::time_point{};
    if(not _trajectorySet)
    {
        throw std::runtime_error(
            "[ERROR] [SERIAL LINK MPCC] step(): No trajectory set. Call set_trajectory() first.");
    }
    if(not std::isfinite(dt) or dt <= 0.0)
    {
        throw std::invalid_argument("[ERROR] [SERIAL LINK MPCC] step(): dt must be positive.");
    }
    _dt = dt;
    update();

    const RobotLibrary::Model::Pose referencePose =
        RobotLibrary::Trajectory::express_pose_in_base(
            _trajectoryFrame, _trajectory.pose_at_progress(_pathProgress));
    const RobotLibrary::Model::Pose currentPose = endpoint_pose();
    const Eigen::Matrix3d referenceRotation = referencePose.quaternion().toRotationMatrix();

    Eigen::Vector<double,ERROR_DIM> error0;
    error0.head<3>() =
        referenceRotation.transpose() * (currentPose.translation() - referencePose.translation());
    const Eigen::Quaterniond relativeOrientation =
        referencePose.quaternion().conjugate() * currentPose.quaternion();
    error0.tail<3>() = quaternion_to_rotation_vector(relativeOrientation);

    _diagnostics.referenceProgress = _pathProgress;
    _diagnostics.error = error0;
    _diagnostics.positionError = error0.head<3>().norm();
    _diagnostics.orientationError = error0.tail<3>().norm();

    const Eigen::Vector<double,NU> optimalControl = solve_mpcc(error0, referenceRotation);
    _uLast = optimalControl;
    _diagnostics.bodyTwist = optimalControl.head<6>();
    _diagnostics.progressRate = optimalControl(6);

    Eigen::Vector<double,6> baseTwist;
    baseTwist.head<3>() = referenceRotation * optimalControl.head<3>();
    baseTwist.tail<3>() = referenceRotation * optimalControl.segment<3>(3);
    _pathProgress = std::clamp(_pathProgress + dt * optimalControl(6), 0.0, 1.0);
    _diagnostics.nextProgress = _pathProgress;
    const auto resolvedRateStart =
        _timingDiagnosticsEnabled ? Clock::now() : Clock::time_point{};
    const Eigen::VectorXd jointCommand = resolve_endpoint_twist(baseTwist);
    if(_timingDiagnosticsEnabled)
    {
        _diagnostics.resolvedRateTimeSeconds =
            std::chrono::duration<double>(Clock::now() - resolvedRateStart).count();
    }
    const Eigen::Vector<double,6> realizedBaseTwist = _jacobianMatrix * jointCommand;
    _diagnostics.realizedBodyTwist.head<3>() =
        referenceRotation.transpose() * realizedBaseTwist.head<3>();
    _diagnostics.realizedBodyTwist.tail<3>() =
        referenceRotation.transpose() * realizedBaseTwist.tail<3>();
    _diagnostics.twistRealizationError =
        (_diagnostics.realizedBodyTwist - _diagnostics.bodyTwist).norm();
    if(_timingDiagnosticsEnabled)
    {
        _diagnostics.totalStepTimeSeconds =
            std::chrono::duration<double>(Clock::now() - totalStart).count();
    }
    ++_controlStepIndex;
    return jointCommand;
}

Eigen::VectorXd
SerialLinkMPCC::step(const double dt, const double estimatedProgress)
{
    if(not std::isfinite(estimatedProgress))
        throw std::invalid_argument("[ERROR] [SERIAL LINK MPCC] step(): estimatedProgress must be finite.");
    _pathProgress = std::clamp(estimatedProgress, 0.0, 1.0);
    return step(dt);
}

Eigen::Vector<double,SerialLinkMPCC::NU>
SerialLinkMPCC::solve_mpcc(const Eigen::Vector<double,ERROR_DIM> &error0,
                           const Eigen::Matrix3d &referenceRotation)
{
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Clock = std::chrono::steady_clock;

    const auto referencePreparationStart =
        _timingDiagnosticsEnabled ? Clock::now() : Clock::time_point{};

    const int N = static_cast<int>(_horizon);
    const int errorDim = ERROR_DIM * N;
    const int controlDim = NU * N;
    const double remaining = std::max(0.0, 1.0 - _pathProgress);
    const double progressMinimum =
        std::min(_vProgressMin,
                 remaining / std::max(static_cast<double>(N) * _dt, 1e-9));
    const Eigen::VectorXd fixedProgressRates = reference_schedule_rates(
        _pathProgress, N, _dt, 1.0 / _vProgressNominal);

    std::vector<Eigen::Vector<double,ERROR_DIM>> pathTangents(static_cast<size_t>(N));
    std::vector<Eigen::Vector<double,ERROR_DIM>> referenceMotionJacobians(
        static_cast<size_t>(N));
    std::vector<Eigen::Vector<double,ERROR_DIM>> referenceMotionOffsets(
        static_cast<size_t>(N), Eigen::Vector<double,ERROR_DIM>::Zero());
    std::vector<Eigen::Vector<double,ERROR_DIM>> referenceMotions(
        static_cast<size_t>(N));
    const bool parentMotionActive = _parentFrameMotion.has_velocity()
        && _parentFrameMotion.body_twist().squaredNorm() > 0.0;
    const auto parentTransform = [&](const int stage)
    {
        return _parentFrameMotion.predicted_pose(stage, _dt);
    };
    double progress = _pathProgress;
    for(int stage = 0; stage < N; ++stage)
    {
        pathTangents[static_cast<size_t>(stage)] =
            path_tangent_at_progress(progress, referenceRotation);

        double rate = fixedProgressRates(stage);
        if(not _fixedProgressSchedule and _warmStart.size() == controlDim)
        {
            rate = _warmStart(stage * NU + 6);
        }
        referenceMotionJacobians[static_cast<size_t>(stage)] =
            pathTangents[static_cast<size_t>(stage)];
        referenceMotions[static_cast<size_t>(stage)] =
            pathTangents[static_cast<size_t>(stage)] * rate;
        if(parentMotionActive)
        {
            RobotLibrary::Model::Pose pathPoseModel =
                _trajectory.pose_at_progress(progress);
            const Eigen::Matrix4d pathPose = pathPoseModel.as_matrix();
            const Eigen::Vector<double,6> bodyTangent =
                _trajectory.tangent_at_progress(progress);
            const Eigen::Matrix3d stageRotation =
                mpcc_stage_reference_rotation(parentTransform(stage), pathPose);
            const auto mappedMotion = [&](const double candidateRate)
            {
                const Eigen::Matrix4d displacement =
                    legacy_repaired_reference_displacement(
                        pathPose, bodyTangent, candidateRate, _dt,
                        parentTransform(stage), parentTransform(stage + 1));
                return mpcc_express_body_tangent_in_prediction_frame(
                    RobotLibrary::Math::se3_logarithm(displacement) / _dt,
                    stageRotation, referenceRotation);
            };
            constexpr double rateStep = 1e-6;
            const Eigen::Vector<double,ERROR_DIM> motion = mappedMotion(rate);
            const Eigen::Vector<double,ERROR_DIM> derivative =
                (mappedMotion(rate + rateStep) - mappedMotion(rate - rateStep))
                / (2.0 * rateStep);
            referenceMotions[static_cast<size_t>(stage)] = motion;
            referenceMotionJacobians[static_cast<size_t>(stage)] = derivative;
            referenceMotionOffsets[static_cast<size_t>(stage)] =
                motion - derivative * rate;
        }
        progress = std::clamp(progress + _dt * (_fixedProgressSchedule
                                  ? rate
                                  : std::clamp(rate, progressMinimum, _vProgressMax)),
                              0.0,
                              1.0);
    }
    if(_timingDiagnosticsEnabled)
    {
        _diagnostics.referencePreparationTimeSeconds =
            std::chrono::duration<double>(Clock::now() - referencePreparationStart).count();
    }

    // Linearised error prediction:
    // e_j = e_0 + dt * sum_{i<=j}(u_i - tau(s_i) * sdot_i).
    MatrixXd errorResponse = MatrixXd::Zero(errorDim, controlDim);
    VectorXd errorOffset = VectorXd::Zero(errorDim);
    MatrixXd errorWeight = MatrixXd::Zero(errorDim, errorDim);
    const auto predictionStart =
        _timingDiagnosticsEnabled ? Clock::now() : Clock::time_point{};
    for(int stage = 0; stage < N; ++stage)
    {
        const int row = stage * ERROR_DIM;
        errorOffset.segment<ERROR_DIM>(row) = error0;

        for(int input = 0; input <= stage; ++input)
        {
            const int column = input * NU;
            errorResponse.block<ERROR_DIM,ERROR_DIM>(row, column) +=
                _dt * Eigen::Matrix<double,ERROR_DIM,ERROR_DIM>::Identity();
            errorResponse.block<ERROR_DIM,1>(row, column + 6) +=
                -_dt * referenceMotionJacobians[static_cast<size_t>(input)];
            errorOffset.segment<ERROR_DIM>(row) -=
                _dt * referenceMotionOffsets[static_cast<size_t>(input)];
        }
    }
    if(_timingDiagnosticsEnabled)
    {
        _diagnostics.errorPredictionTimeSeconds =
            std::chrono::duration<double>(Clock::now() - predictionStart).count();
    }

    const auto costWeightStart =
        _timingDiagnosticsEnabled ? Clock::now() : Clock::time_point{};
    for(int stage = 0; stage < N; ++stage)
    {
        const int row = stage * ERROR_DIM;
        const Eigen::Vector<double,ERROR_DIM> &tangent =
            pathTangents[static_cast<size_t>(stage)];
        const Eigen::Matrix3d positionLagProjection =
            rmpcc_decoupled_error_projection(tangent, 1e-9).lag.block<3,3>(0,0);

        double contourWeight = _wContour;
        double lagWeight = _wLag;
        double orientationWeight = _wOrientation;
        if(_ablationProfile == MpccAblationProfile::MatchedFeedbackGain
           or _ablationProfile == MpccAblationProfile::OptimizedProgress)
        {
            contourWeight = RobotLibrary::Math::integrator_stage_weight_for_gain(
                _matchedPositionGain, 2.5, _dt);
            orientationWeight = RobotLibrary::Math::integrator_stage_weight_for_gain(
                _matchedOrientationGain, 1.6, _dt);
            if(stage == N - 1)
            {
                contourWeight = RobotLibrary::Math::integrator_terminal_weight_for_gain(
                    _matchedPositionGain, 2.5, _dt);
                lagWeight = RobotLibrary::Math::integrator_terminal_weight(
                    _wLag, 2.5, _dt);
                orientationWeight = RobotLibrary::Math::integrator_terminal_weight_for_gain(
                    _matchedOrientationGain, 1.6, _dt);
            }
        }
        const Eigen::Matrix3d positionWeight =
            contourWeight * Eigen::Matrix3d::Identity()
            + (lagWeight - contourWeight) * positionLagProjection;
        errorWeight.block<3,3>(row, row) = positionWeight;
        errorWeight.block<3,3>(row + 3, row + 3) =
            orientationWeight * Eigen::Matrix3d::Identity();
    }
    if(_timingDiagnosticsEnabled)
    {
        _diagnostics.costWeightTimeSeconds =
            std::chrono::duration<double>(Clock::now() - costWeightStart).count();
    }

    const auto pathVelocityStart =
        _timingDiagnosticsEnabled ? Clock::now() : Clock::time_point{};
    const bool feedforwardCorrection =
        _ablationProfile == MpccAblationProfile::FeedforwardCorrection;
    const bool optimizedProgress =
        _ablationProfile == MpccAblationProfile::OptimizedProgress;
    const bool matchedFeedbackGain =
        _ablationProfile == MpccAblationProfile::MatchedFeedbackGain;
    // The fixed-progress profile uses exact discrete-LQR stage and terminal
    // weights above.  This avoids the horizon-dependent sqrt(Q/R) shortcut.
    const double inputLinearWeight =
        (feedforwardCorrection or optimizedProgress or matchedFeedbackGain)
        ? 0.0 : _wInputLinear;
    const double inputAngularWeight =
        (feedforwardCorrection or optimizedProgress or matchedFeedbackGain)
        ? 0.0 : _wInputAngular;
    const double inputProgressWeight =
        (feedforwardCorrection or optimizedProgress) ? 0.0 : _wInputProgress;
    Eigen::Matrix<double,ERROR_DIM,ERROR_DIM> pathVelocityWeight =
        _wVelocityTracking * Eigen::Matrix<double,ERROR_DIM,ERROR_DIM>::Identity();
    if(feedforwardCorrection)
    {
        // Penalise the correction u - tau(s) sdot.  Penalising u and sdot
        // separately makes the optimum prefer a stopped trajectory even when
        // the geometric tracking error is zero.
        pathVelocityWeight.diagonal().head<3>().array() += _wInputLinear;
        pathVelocityWeight.diagonal().tail<3>().array() += _wInputAngular;
    }
    else if(matchedFeedbackGain or optimizedProgress)
    {
        // Preserve path feedforward and give the correction the same local
        // bandwidth in fixed- and optimized-progress experiments.
        pathVelocityWeight.diagonal() << 2.5, 2.5, 2.5, 1.6, 1.6, 1.6;
    }
    const MatrixXd pathVelocityWeightHorizon =
        RobotLibrary::Math::block_diagonal(pathVelocityWeight, static_cast<unsigned int>(N));

    Eigen::Matrix<double,NU,NU> inputWeight =
        Eigen::Matrix<double,NU,NU>::Zero();
    inputWeight.diagonal() <<
        inputLinearWeight, inputLinearWeight, inputLinearWeight,
        inputAngularWeight, inputAngularWeight, inputAngularWeight,
        inputProgressWeight;
    const MatrixXd inputWeightHorizon =
        RobotLibrary::Math::block_diagonal(inputWeight, static_cast<unsigned int>(N));

    VectorXd nominalControl = VectorXd::Zero(controlDim);
    MatrixXd pathVelocityResponse = MatrixXd::Zero(errorDim, controlDim);
    VectorXd pathVelocityOffset = VectorXd::Zero(errorDim);
    const double nominalSeedRate =
        std::min(_vProgressNominal,
                 remaining / std::max(static_cast<double>(N) * _dt, 1e-9));
    for(int stage = 0; stage < N; ++stage)
    {
        const int offset = stage * NU;
        const int row = stage * ERROR_DIM;
        const double stageSeedRate = _fixedProgressSchedule
            ? fixedProgressRates(stage) : nominalSeedRate;
        nominalControl.segment<ERROR_DIM>(offset) =
            parentMotionActive
            ? referenceMotions[static_cast<size_t>(stage)]
            : pathTangents[static_cast<size_t>(stage)] * stageSeedRate;
        nominalControl(offset + 6) = stageSeedRate;
        pathVelocityResponse.block<ERROR_DIM,ERROR_DIM>(row, offset).setIdentity();
        pathVelocityResponse.block<ERROR_DIM,1>(row, offset + 6) =
            -referenceMotionJacobians[static_cast<size_t>(stage)];
        pathVelocityOffset.segment<ERROR_DIM>(row) =
            -referenceMotionOffsets[static_cast<size_t>(stage)];
    }
    if(_timingDiagnosticsEnabled)
    {
        _diagnostics.pathVelocityObjectiveTimeSeconds =
            std::chrono::duration<double>(Clock::now() - pathVelocityStart).count();
    }

    const auto hessianStart =
        _timingDiagnosticsEnabled ? Clock::now() : Clock::time_point{};
    MatrixXd difference = MatrixXd::Zero(controlDim, controlDim);
    for(int stage = 0; stage < N; ++stage)
    {
        difference.block<NU,NU>(stage * NU, stage * NU).setIdentity();
        if(stage > 0)
        {
            difference.block<NU,NU>(stage * NU, (stage - 1) * NU) =
                -Eigen::Matrix<double,NU,NU>::Identity();
        }
    }
    const MatrixXd deltaWeight =
        _wDeltaU * MatrixXd::Identity(controlDim, controlDim);
    VectorXd deltaReference = VectorXd::Zero(controlDim);
    deltaReference.head<NU>() = _uLast;

    MatrixXd H =
        errorResponse.transpose() * errorWeight * errorResponse
        + inputWeightHorizon
        + difference.transpose() * deltaWeight * difference
        + pathVelocityResponse.transpose() * pathVelocityWeightHorizon * pathVelocityResponse;
    H += 1e-8 * MatrixXd::Identity(controlDim, controlDim);
    H = 0.5 * (H + H.transpose());

    VectorXd f =
        errorResponse.transpose() * errorWeight * errorOffset
        + pathVelocityResponse.transpose() * pathVelocityWeightHorizon
            * pathVelocityOffset
        - difference.transpose() * deltaWeight * deltaReference;
    const double progressReward = optimizedProgress ? 0.2 : _qProgressReward;
    for(int stage = 0; stage < N; ++stage)
    {
        f(stage * NU + 6) -= progressReward * _dt;
    }
    if(_timingDiagnosticsEnabled)
    {
        _diagnostics.hessianAssemblyTimeSeconds =
            std::chrono::duration<double>(Clock::now() - hessianStart).count();
    }

    const auto constraintStart =
        _timingDiagnosticsEnabled ? Clock::now() : Clock::time_point{};
    VectorXd lower = VectorXd::Zero(controlDim);
    VectorXd upper = VectorXd::Zero(controlDim);
    for(int stage = 0; stage < N; ++stage)
    {
        const int offset = stage * NU;
        lower.segment<3>(offset).setConstant(-_vMaxLinear);
        upper.segment<3>(offset).setConstant(_vMaxLinear);
        lower.segment<3>(offset + 3).setConstant(-_vMaxAngular);
        upper.segment<3>(offset + 3).setConstant(_vMaxAngular);
        lower(offset + 6) = _fixedProgressSchedule
            ? fixedProgressRates(stage) : progressMinimum;
        upper(offset + 6) = _fixedProgressSchedule
            ? fixedProgressRates(stage) : _vProgressMax;
    }

    const RobotLibrary::Math::BoxConstraint box =
        RobotLibrary::Math::box_constraint(lower, upper);
    MatrixXd constraintMatrix =
        MatrixXd::Zero(box.constraintMatrix.rows() + 1, controlDim);
    VectorXd constraintVector =
        VectorXd::Zero(box.constraintVector.size() + 1);
    constraintMatrix.topRows(box.constraintMatrix.rows()) = box.constraintMatrix;
    constraintVector.head(box.constraintVector.size()) = box.constraintVector;
    for(int stage = 0; stage < N; ++stage)
    {
        constraintMatrix(box.constraintMatrix.rows(), stage * NU + 6) = _dt;
    }
    constraintVector.tail<1>()(0) = remaining;

    auto make_feasible = [&](VectorXd seed)
    {
        if(seed.size() != controlDim)
        {
            seed = nominalControl;
        }
        for(int i = 0; i < controlDim; ++i)
        {
            seed(i) = std::clamp(seed(i), lower(i), upper(i));
        }

        double predictedAdvance = 0.0;
        for(int stage = 0; stage < N; ++stage)
        {
            predictedAdvance += _dt * seed(stage * NU + 6);
        }
        if(not _fixedProgressSchedule and predictedAdvance > remaining)
        {
            const double feasibleRate =
                (N > 0 and _dt > 0.0) ? remaining / (static_cast<double>(N) * _dt) : 0.0;
            for(int stage = 0; stage < N; ++stage)
            {
                seed(stage * NU + 6) = feasibleRate;
            }
        }
        return seed;
    };

    VectorXd seed = make_feasible(_warmStart);
    if(_timingDiagnosticsEnabled)
    {
        _diagnostics.constraintConstructionTimeSeconds =
            std::chrono::duration<double>(Clock::now() - constraintStart).count();
    }
    _diagnostics.qpStatus = 0.0;
    const auto qpStart = _timingDiagnosticsEnabled ? Clock::now() : Clock::time_point{};
    VectorXd optimum;
    SolverResults<double> qpResults;
    if(_fixedProgressSchedule)
    {
        optimum = _qpSolver.solve(H, f, constraintMatrix, constraintVector, seed);
        qpResults = _qpSolver.results();
    }
    else
    {
        const auto specialized = solve_box_aware_active_set(
            H, f, constraintMatrix, constraintVector, seed, _qpOptions);
        optimum = specialized.solution;
        qpResults = specialized.solver;
    }
    if(_timingDiagnosticsEnabled)
    {
        _diagnostics.qpSolveTimeSeconds =
            std::chrono::duration<double>(Clock::now() - qpStart).count();
    }
    const auto postQpStart =
        _timingDiagnosticsEnabled ? Clock::now() : Clock::time_point{};
    const MatrixXd emptyA(0, controlDim);
    const VectorXd emptyY(0);
    write_realtime_recovery_qp_capture(
        "mpcc", _controlStepIndex, H, f, emptyA, emptyY,
        constraintMatrix, constraintVector, seed, optimum, qpResults);
    _diagnostics.qpIterations = static_cast<double>(qpResults.numberOfSteps);
    _diagnostics.qpFinalStepSize = qpResults.finalStepSize;
    _diagnostics.qpObjective = qpResults.objectiveFunction;
    _diagnostics.qpConverged =
        qpResults.terminationReason == SolverTerminationReason::Converged;
    _diagnostics.qpHitMaxIterations =
        qpResults.terminationReason == SolverTerminationReason::MaxIterations;
    _diagnostics.qpActiveSetChanges = static_cast<double>(qpResults.activeSetChanges);
    _diagnostics.qpMaximumActiveSetSize =
        static_cast<double>(qpResults.maximumActiveSetSize);
    _diagnostics.qpUniqueActiveConstraints =
        static_cast<double>(qpResults.uniqueActiveConstraints);
    if(optimum.size() != controlDim or not optimum.allFinite())
    {
        throw std::runtime_error(
            "[ERROR] [SERIAL LINK MPCC] solve_mpcc(): QP returned an invalid solution.");
    }
    const double constraintViolation =
        (constraintMatrix * optimum - constraintVector).maxCoeff();
    if(not std::isfinite(constraintViolation) or constraintViolation > 1e-6)
    {
        throw std::runtime_error(
            "[ERROR] [SERIAL LINK MPCC] solve_mpcc(): QP returned an infeasible solution (maximum violation "
            + std::to_string(constraintViolation) + ").");
    }
    _diagnostics.qpStatus = 1.0;

    _warmStart = optimum;
    for(int stage = 0; stage < N - 1; ++stage)
    {
        _warmStart.segment<NU>(stage * NU) = optimum.segment<NU>((stage + 1) * NU);
    }

    Eigen::Vector<double,NU> firstControl = optimum.head<NU>();
    for(int i = 0; i < 3; ++i)
    {
        firstControl(i) = std::clamp(firstControl(i), -_vMaxLinear, _vMaxLinear);
        firstControl(i + 3) =
            std::clamp(firstControl(i + 3), -_vMaxAngular, _vMaxAngular);
    }
    firstControl(6) = std::clamp(firstControl(6), lower(6), upper(6));
    _diagnostics.parentFrameMotionActive = parentMotionActive;
    _diagnostics.parentFrameBodyTwist = _parentFrameMotion.body_twist();
    _diagnostics.measuredParentPose = _parentFrameMotion.current_pose();
    _diagnostics.parentMeasurementTimeSeconds = _parentFrameMotion.current_time();
    _diagnostics.predictedParentPoseFirst = parentTransform(1);
    _diagnostics.predictedParentPoseHorizon = parentTransform(N);
    const Eigen::Matrix4d firstPathPose =
        _trajectory.pose_at_progress(_pathProgress).as_matrix();
    _diagnostics.parentReferenceFactorFirst = parent_frame_reference_factor(
        firstPathPose, parentTransform(0), parentTransform(1));
    _diagnostics.repairedReferenceDisplacementFirst =
        legacy_repaired_reference_displacement(
            firstPathPose,
            _trajectory.tangent_at_progress(_pathProgress),
            firstControl(6), _dt, parentTransform(0), parentTransform(1));
    if(_timingDiagnosticsEnabled)
    {
        _diagnostics.postQpTimeSeconds =
            std::chrono::duration<double>(Clock::now() - postQpStart).count();
    }
    return firstControl;
}

Eigen::Vector<double,6>
SerialLinkMPCC::path_tangent_at_progress(
    const double progress,
    const Eigen::Matrix3d &referenceRotation)
{
    // A rigid left multiplication by the trajectory frame leaves the body
    // tangent invariant, but every horizon stage must be expressed in the
    // single prediction frame used by error0 and all task-space controls.
    const Eigen::Matrix3d frameRotation =
        _trajectoryFrame.transformInBase.block<3,3>(0,0);
    const Eigen::Matrix3d stageRotation =
        frameRotation
        * _trajectory.pose_at_progress(progress).quaternion().toRotationMatrix();
    return mpcc_express_body_tangent_in_prediction_frame(
        _trajectory.tangent_at_progress(progress), stageRotation, referenceRotation);
}

} } // namespace
