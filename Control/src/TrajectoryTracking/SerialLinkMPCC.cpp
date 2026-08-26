/**
 * @file    SerialLinkMPCC.cpp
 * @brief   MPCC controller implementation for serial link robot arms.
 */

#include <Control/TrajectoryTracking/SerialLinkMPCC.h>
#include "detail/ProgressSchedule.h"
#include "detail/RmpccCostGeometry.h"
#include <Math/CondensedMPC.h>
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
  _qpSolver(parameters.mpcc_qp_options()),
  _qpStepSizeTolerance(parameters.mpcc_qp_options().stepSizeTolerance)
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
    reset_additional_virtual_state();
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
SerialLinkMPCC::set_objective_parameters(const MpccObjectiveParameters &parameters)
{
    const double values[] = {
        parameters.contourWeight, parameters.lagWeight,
        parameters.orientationWeight, parameters.inputLinearWeight,
        parameters.inputAngularWeight, parameters.inputProgressWeight,
        parameters.inputDifferenceWeight, parameters.pathVelocityWeight,
        parameters.progressReward};
    for(const double value : values)
    {
        if(not std::isfinite(value) or value < 0.0)
        {
            throw std::invalid_argument(
                "[ERROR] [SERIAL LINK MPCC] Objective weights must be finite and nonnegative.");
        }
    }
    _wContour = parameters.contourWeight;
    _wLag = parameters.lagWeight;
    _wOrientation = parameters.orientationWeight;
    _wInputLinear = parameters.inputLinearWeight;
    _wInputAngular = parameters.inputAngularWeight;
    _wInputProgress = parameters.inputProgressWeight;
    _wDeltaU = parameters.inputDifferenceWeight;
    _wVelocityTracking = parameters.pathVelocityWeight;
    _qProgressReward = parameters.progressReward;
}

void
SerialLinkMPCC::set_progress_rate_limits(
    const double minimum,
    const double nominal,
    const double maximum)
{
    if(not std::isfinite(minimum) or not std::isfinite(nominal)
       or not std::isfinite(maximum) or minimum < 0.0
       or minimum > nominal or nominal > maximum)
    {
        throw std::invalid_argument(
            "[ERROR] [SERIAL LINK MPCC] Progress rates must satisfy 0 <= min <= nominal <= max.");
    }
    _vProgressMin = minimum;
    _vProgressNominal = nominal;
    _vProgressMax = maximum;
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
    return step_impl(dt, 0.0);
}

Eigen::VectorXd
SerialLinkMPCC::step_at_time(const double controlTimeSeconds, const double dt)
{
    if(not std::isfinite(controlTimeSeconds))
    {
        throw std::invalid_argument(
            "[ERROR] [SERIAL LINK MPCC] step_at_time(): control time must be finite.");
    }
    const double measurementAge = std::max(
        0.0, controlTimeSeconds - _parentFrameMotion.current_time());
    return step_impl(dt, measurementAge);
}

Eigen::VectorXd
SerialLinkMPCC::step_at_time(
    const double controlTimeSeconds,
    const double dt,
    const double estimatedProgress)
{
    if(not std::isfinite(estimatedProgress))
    {
        throw std::invalid_argument(
            "[ERROR] [SERIAL LINK MPCC] step_at_time(): estimated progress must be finite.");
    }
    _pathProgress = std::clamp(estimatedProgress, 0.0, 1.0);
    return step_at_time(controlTimeSeconds, dt);
}

Eigen::VectorXd
SerialLinkMPCC::step_impl(const double dt, const double parentMeasurementAgeSeconds)
{
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
        reference_pose_at_progress(_pathProgress);
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

    const Eigen::Vector<double,BASE_STAGE_CONTROL_DIMENSION> optimalControl = solve_mpcc(
        error0, referenceRotation, currentPose.translation(),
        currentPose.quaternion().toRotationMatrix(), parentMeasurementAgeSeconds);
    _uLast = optimalControl;
    _diagnostics.bodyTwist = optimalControl.head<6>();
    _diagnostics.progressRate = optimalControl(6);

    Eigen::Vector<double,6> baseTwist;
    baseTwist.head<3>() = referenceRotation * optimalControl.head<3>();
    baseTwist.tail<3>() = referenceRotation * optimalControl.segment<3>(3);
    baseTwist = postprocess_base_twist(baseTwist, dt);
    if(not baseTwist.allFinite())
    {
        throw std::runtime_error(
            "[ERROR] [SERIAL LINK MPCC] step(): Cartesian postprocessor returned a non-finite twist.");
    }
    _diagnostics.commandedBaseTwist = baseTwist;
    _pathProgress = std::clamp(_pathProgress + dt * optimalControl(6), 0.0, 1.0);
    _diagnostics.nextProgress = _pathProgress;
    const Eigen::VectorXd jointCommand = resolve_endpoint_twist(baseTwist);
    const Eigen::Vector<double,6> realizedBaseTwist = _jacobianMatrix * jointCommand;
    _diagnostics.realizedBaseTwist = realizedBaseTwist;
    _diagnostics.realizedBodyTwist.head<3>() =
        referenceRotation.transpose() * realizedBaseTwist.head<3>();
    _diagnostics.realizedBodyTwist.tail<3>() =
        referenceRotation.transpose() * realizedBaseTwist.tail<3>();
    _diagnostics.twistRealizationError =
        (realizedBaseTwist - baseTwist).norm();
    on_twist_resolved(baseTwist, realizedBaseTwist, dt);
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

Eigen::Vector<double,SerialLinkMPCC::BASE_STAGE_CONTROL_DIMENSION>
SerialLinkMPCC::solve_mpcc(const Eigen::Vector<double,ERROR_DIM> &error0,
                           const Eigen::Matrix3d &referenceRotation,
                           const Eigen::Vector3d &currentToolPositionBase,
                           const Eigen::Matrix3d &currentToolRotationBase,
                           const double parentMeasurementAgeSeconds)
{
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    const int N = static_cast<int>(_horizon);
    const int errorDim = ERROR_DIM * N;
    const int stageControlDim = stage_control_dimension();
    if(stageControlDim < BASE_STAGE_CONTROL_DIMENSION)
    {
        throw std::logic_error(
            "[ERROR] [SERIAL LINK MPCC] stage_control_dimension() must be at least 7.");
    }
    const int controlDim = stageControlDim * N;
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
    std::vector<Eigen::MatrixXd> additionalReferenceTangents(
        static_cast<size_t>(N));
    std::vector<double> stageProgress(static_cast<size_t>(N), _pathProgress);
    const bool parentMotionActive = _parentFrameMotion.has_velocity()
        && _parentFrameMotion.body_twist().squaredNorm() > 0.0;
    std::vector<Eigen::Matrix4d> parentTransforms(static_cast<size_t>(N + 1));
    for(int stage = 0; stage <= N; ++stage)
    {
        parentTransforms[static_cast<size_t>(stage)] =
            _parentFrameMotion.predicted_pose_with_age(
                parentMeasurementAgeSeconds, stage, _dt);
    }
    const auto parentTransform = [&](const int stage) -> const Eigen::Matrix4d &
    {
        return parentTransforms[static_cast<size_t>(stage)];
    };
    double progress = _pathProgress;
    for(int stage = 0; stage < N; ++stage)
    {
        stageProgress[static_cast<size_t>(stage)] = progress;
        pathTangents[static_cast<size_t>(stage)] =
            path_tangent_at_progress(progress, referenceRotation);

        additionalReferenceTangents[static_cast<size_t>(stage)] =
            additional_reference_tangents(
                stage, progress, referenceRotation, parentTransform(stage));
        if(additionalReferenceTangents[static_cast<size_t>(stage)].rows() != ERROR_DIM
           or additionalReferenceTangents[static_cast<size_t>(stage)].cols()
                != stageControlDim - BASE_STAGE_CONTROL_DIMENSION
           or not additionalReferenceTangents[static_cast<size_t>(stage)].allFinite())
        {
            throw std::runtime_error(
                "[ERROR] [SERIAL LINK MPCC] Invalid additional reference tangent dimensions.");
        }

        double rate = fixedProgressRates(stage);
        if(not _fixedProgressSchedule and _warmStart.size() >= controlDim)
        {
            rate = _warmStart(stage * stageControlDim + PROGRESS_RATE_INDEX);
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

    // Linearised error prediction:
    // e_j = e_0 + dt * sum_{i<=j}(xi_i - tau_s(s_i)*sdot_i
    //                              - tau_virtual(s_i)*virtual_rate_i).
    MatrixXd errorResponse = MatrixXd::Zero(errorDim, controlDim);
    VectorXd errorOffset = VectorXd::Zero(errorDim);
    MatrixXd errorWeight = MatrixXd::Zero(errorDim, errorDim);
    for(int stage = 0; stage < N; ++stage)
    {
        const int row = stage * ERROR_DIM;
        errorOffset.segment<ERROR_DIM>(row) = error0;

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
        const Eigen::Matrix3d defaultPositionWeight =
            contourWeight * Eigen::Matrix3d::Identity()
            + (lagWeight - contourWeight) * positionLagProjection;
        const Eigen::Matrix3d positionWeight = position_error_weight(
            stage, tangent, defaultPositionWeight, referenceRotation,
            parentTransform(stage + 1));
        if(not positionWeight.allFinite())
        {
            throw std::runtime_error(
                "[ERROR] [SERIAL LINK MPCC] position_error_weight(): non-finite matrix.");
        }
        errorWeight.block<3,3>(row, row) = positionWeight;
        errorWeight.block<3,3>(row + 3, row + 3) =
            orientationWeight * Eigen::Matrix3d::Identity();

        for(int input = 0; input <= stage; ++input)
        {
            const int column = input * stageControlDim;
            errorResponse.block<ERROR_DIM,ERROR_DIM>(row, column) +=
                _dt * Eigen::Matrix<double,ERROR_DIM,ERROR_DIM>::Identity();
            errorResponse.block<ERROR_DIM,1>(row, column + PROGRESS_RATE_INDEX) +=
                -_dt * referenceMotionJacobians[static_cast<size_t>(input)];
            const int additionalCount =
                stageControlDim - BASE_STAGE_CONTROL_DIMENSION;
            if(additionalCount > 0)
            {
                errorResponse.block(
                    row, column + BASE_STAGE_CONTROL_DIMENSION,
                    ERROR_DIM, additionalCount) +=
                    -_dt * additionalReferenceTangents[static_cast<size_t>(input)];
            }
            errorOffset.segment<ERROR_DIM>(row) -=
                _dt * referenceMotionOffsets[static_cast<size_t>(input)];
        }
    }

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

    MatrixXd inputWeight = MatrixXd::Zero(stageControlDim, stageControlDim);
    inputWeight.diagonal().head<BASE_STAGE_CONTROL_DIMENSION>() <<
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
        const int offset = stage * stageControlDim;
        const int row = stage * ERROR_DIM;
        const double stageSeedRate = _fixedProgressSchedule
            ? fixedProgressRates(stage) : nominalSeedRate;
        nominalControl.segment<ERROR_DIM>(offset) =
            parentMotionActive
            ? referenceMotions[static_cast<size_t>(stage)]
            : pathTangents[static_cast<size_t>(stage)] * stageSeedRate;
        nominalControl(offset + PROGRESS_RATE_INDEX) = stageSeedRate;
        pathVelocityResponse.block<ERROR_DIM,ERROR_DIM>(row, offset).setIdentity();
        pathVelocityResponse.block<ERROR_DIM,1>(row, offset + PROGRESS_RATE_INDEX) =
            -referenceMotionJacobians[static_cast<size_t>(stage)];
        const int additionalCount =
            stageControlDim - BASE_STAGE_CONTROL_DIMENSION;
        if(additionalCount > 0)
        {
            pathVelocityResponse.block(
                row, offset + BASE_STAGE_CONTROL_DIMENSION,
                ERROR_DIM, additionalCount) =
                -additionalReferenceTangents[static_cast<size_t>(stage)];
        }
        pathVelocityOffset.segment<ERROR_DIM>(row) =
            -referenceMotionOffsets[static_cast<size_t>(stage)];
    }

    const int baseDifferenceDim = BASE_STAGE_CONTROL_DIMENSION * N;
    MatrixXd difference = MatrixXd::Zero(baseDifferenceDim, controlDim);
    for(int stage = 0; stage < N; ++stage)
    {
        difference.block(
            stage * BASE_STAGE_CONTROL_DIMENSION,
            stage * stageControlDim,
            BASE_STAGE_CONTROL_DIMENSION,
            BASE_STAGE_CONTROL_DIMENSION).setIdentity();
        if(stage > 0)
        {
            difference.block(
                stage * BASE_STAGE_CONTROL_DIMENSION,
                (stage - 1) * stageControlDim,
                BASE_STAGE_CONTROL_DIMENSION,
                BASE_STAGE_CONTROL_DIMENSION) =
                -Eigen::Matrix<double,BASE_STAGE_CONTROL_DIMENSION,
                               BASE_STAGE_CONTROL_DIMENSION>::Identity();
        }
    }
    const MatrixXd deltaWeight =
        _wDeltaU * MatrixXd::Identity(baseDifferenceDim, baseDifferenceDim);
    VectorXd deltaReference = VectorXd::Zero(baseDifferenceDim);
    deltaReference.head<BASE_STAGE_CONTROL_DIMENSION>() = _uLast;

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
        f(stage * stageControlDim + PROGRESS_RATE_INDEX) -= progressReward * _dt;
    }

    VectorXd lower = VectorXd::Zero(controlDim);
    VectorXd upper = VectorXd::Zero(controlDim);
    for(int stage = 0; stage < N; ++stage)
    {
        const int offset = stage * stageControlDim;
        lower.segment<3>(offset).setConstant(-_vMaxLinear);
        upper.segment<3>(offset).setConstant(_vMaxLinear);
        lower.segment<3>(offset + 3).setConstant(-_vMaxAngular);
        upper.segment<3>(offset + 3).setConstant(_vMaxAngular);
        lower(offset + PROGRESS_RATE_INDEX) = _fixedProgressSchedule
            ? fixedProgressRates(stage) : progressMinimum;
        upper(offset + PROGRESS_RATE_INDEX) = _fixedProgressSchedule
            ? fixedProgressRates(stage) : _vProgressMax;
        configure_additional_stage_inputs(stage, lower, upper, nominalControl);
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
        constraintMatrix(box.constraintMatrix.rows(),
                         stage * stageControlDim + PROGRESS_RATE_INDEX) = _dt;
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
            predictedAdvance += _dt * seed(
                stage * stageControlDim + PROGRESS_RATE_INDEX);
        }
        if(not _fixedProgressSchedule and predictedAdvance > remaining)
        {
            const double feasibleRate =
                (N > 0 and _dt > 0.0) ? remaining / (static_cast<double>(N) * _dt) : 0.0;
            for(int stage = 0; stage < N; ++stage)
            {
                seed(stage * stageControlDim + PROGRESS_RATE_INDEX) = feasibleRate;
            }
        }
        return seed;
    };

    VectorXd baseWarmStart;
    if(_warmStart.size() >= controlDim)
    {
        baseWarmStart = _warmStart.head(controlDim);
    }
    VectorXd seed = make_feasible(baseWarmStart);
    MpccQpExtensionContext extensionContext;
    extensionContext.horizon = N;
    extensionContext.baseControlDimension = controlDim;
    extensionContext.stageControlDimension = stageControlDim;
    extensionContext.dt = _dt;
    extensionContext.predictionRotation = referenceRotation;
    extensionContext.currentToolPositionBase = currentToolPositionBase;
    extensionContext.currentToolRotationBase = currentToolRotationBase;
    extensionContext.parentTransforms = parentTransforms;
    extensionContext.stageProgress = stageProgress;
    extensionContext.previousWarmStart = _warmStart;
    extend_qp_problem(extensionContext, H, f, constraintMatrix, constraintVector, seed);
    if(H.rows() != H.cols() or f.size() != H.rows()
       or constraintMatrix.cols() != H.cols()
       or constraintVector.size() != constraintMatrix.rows()
       or seed.size() != H.rows() or not H.allFinite() or not f.allFinite())
    {
        throw std::runtime_error(
            "[ERROR] [SERIAL LINK MPCC] extend_qp_problem(): inconsistent augmented QP.");
    }
    H = 0.5 * (H + H.transpose());
    _diagnostics.qpHessian = H;
    _diagnostics.qpGradient = f;
    _diagnostics.qpConstraintMatrix = constraintMatrix;
    _diagnostics.qpConstraintVector = constraintVector;
    _diagnostics.qpSeed = seed;
    _diagnostics.qpStatus = 0.0;
    const auto solveStart = std::chrono::steady_clock::now();
    VectorXd optimum = _qpSolver.solve(H, f, constraintMatrix, constraintVector, seed);
    _diagnostics.qpSolveTimeSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - solveStart).count();
    const SolverResults<double> qpResults = _qpSolver.results();
    _diagnostics.qpConverged = qpResults.converged;
    _diagnostics.qpIterations = static_cast<double>(qpResults.numberOfSteps);
    _diagnostics.qpFinalStepSize = qpResults.finalStepSize;
    _diagnostics.qpObjective = qpResults.objectiveFunction;
    _diagnostics.qpReturnedSolution = optimum;
    if(optimum.size() != H.rows() or not optimum.allFinite())
    {
        throw std::runtime_error(
            "[ERROR] [SERIAL LINK MPCC] solve_mpcc(): QP returned an invalid solution.");
    }
    double constraintViolation =
        (constraintMatrix * optimum - constraintVector).maxCoeff();
    _diagnostics.qpMaximumConstraintViolation = constraintViolation;
    // The active-set stopping test is intentionally frozen at 1e-4 for the
    // outer MPCC.  It can therefore terminate at the correct active face with
    // an O(1e-6) floating-point half-space residual.  Restore feasibility only
    // for residuals already inside that solver tolerance; larger violations
    // remain hard failures and are never hidden by this numerical polish.
    if(std::isfinite(constraintViolation)
       && constraintViolation > 1e-9
       && constraintViolation <= _qpStepSizeTolerance)
    {
        for(int sweep = 0; sweep < 20 && constraintViolation > 1e-9; ++sweep)
        {
            for(int row = 0; row < constraintMatrix.rows(); ++row)
            {
                const double residual =
                    constraintMatrix.row(row).dot(optimum) - constraintVector(row);
                const double normSquared = constraintMatrix.row(row).squaredNorm();
                if(residual > 1e-10 && normSquared > 1e-20)
                {
                    optimum -= constraintMatrix.row(row).transpose()
                        * ((residual + 1e-10) / normSquared);
                }
            }
            constraintViolation =
                (constraintMatrix * optimum - constraintVector).maxCoeff();
            _diagnostics.qpMaximumConstraintViolation = constraintViolation;
        }
        _diagnostics.qpObjective = 0.5 * optimum.dot(H * optimum) + f.dot(optimum);
    }
    if(not std::isfinite(constraintViolation) or constraintViolation > 1e-6)
    {
        const Eigen::VectorXd seedResidual =
            constraintMatrix * seed - constraintVector;
        Eigen::Index seedViolationRow = -1;
        const double seedViolation = seedResidual.maxCoeff(&seedViolationRow);
        throw std::runtime_error(
            "[ERROR] [SERIAL LINK MPCC] solve_mpcc(): QP returned an infeasible solution (maximum violation "
            + std::to_string(constraintViolation)
            + ", seed violation " + std::to_string(seedViolation)
            + " at row " + std::to_string(seedViolationRow)
            + " (lhs "
            + std::to_string(constraintMatrix.row(seedViolationRow).dot(seed))
            + ", rhs " + std::to_string(constraintVector(seedViolationRow)) + ")"
            + ", iterations " + std::to_string(qpResults.numberOfSteps)
            + ", final step " + std::to_string(qpResults.finalStepSize)
            + ", solution norm " + std::to_string(optimum.norm()) + ").");
    }
    _diagnostics.qpStatus = 1.0;
    _diagnostics.stageControlDimension = stageControlDim;
    _diagnostics.optimalHorizon = optimum;

    _warmStart = optimum;
    for(int stage = 0; stage < N - 1; ++stage)
    {
        _warmStart.segment(stage * stageControlDim, stageControlDim) =
            optimum.segment((stage + 1) * stageControlDim, stageControlDim);
    }
    shift_extension_warm_start(extensionContext, optimum, _warmStart);
    _diagnostics.shiftedWarmStart = _warmStart;
    on_extended_qp_solution(extensionContext, optimum);

    Eigen::Vector<double,BASE_STAGE_CONTROL_DIMENSION> firstControl =
        optimum.head<BASE_STAGE_CONTROL_DIMENSION>();
    for(int i = 0; i < 3; ++i)
    {
        firstControl(i) = std::clamp(firstControl(i), -_vMaxLinear, _vMaxLinear);
        firstControl(i + 3) =
            std::clamp(firstControl(i + 3), -_vMaxAngular, _vMaxAngular);
    }
    firstControl(PROGRESS_RATE_INDEX) = std::clamp(
        firstControl(PROGRESS_RATE_INDEX),
        lower(PROGRESS_RATE_INDEX), upper(PROGRESS_RATE_INDEX));
    _diagnostics.parentFrameMotionActive = parentMotionActive;
    _diagnostics.parentFrameBodyTwist = _parentFrameMotion.body_twist();
    _diagnostics.measuredParentPose = _parentFrameMotion.current_pose();
    _diagnostics.parentMeasurementTimeSeconds = _parentFrameMotion.current_time();
    _diagnostics.parentMeasurementAgeSeconds = parentMeasurementAgeSeconds;
    _diagnostics.predictedParentPoseFirst = parentTransform(1);
    _diagnostics.predictedParentPoseHorizon = parentTransform(N);
    _diagnostics.predictedParentTransforms = parentTransforms;
    const Eigen::Matrix4d firstPathPose =
        _trajectory.pose_at_progress(_pathProgress).as_matrix();
    _diagnostics.parentReferenceFactorFirst = parent_frame_reference_factor(
        firstPathPose, parentTransform(0), parentTransform(1));
    _diagnostics.repairedReferenceDisplacementFirst =
        legacy_repaired_reference_displacement(
            firstPathPose,
            _trajectory.tangent_at_progress(_pathProgress),
            firstControl(6), _dt, parentTransform(0), parentTransform(1));
    return firstControl;
}

RobotLibrary::Model::Pose
SerialLinkMPCC::reference_pose_at_progress(const double progress)
{
    return RobotLibrary::Trajectory::express_pose_in_base(
        _trajectoryFrame, _trajectory.pose_at_progress(progress));
}

Eigen::MatrixXd
SerialLinkMPCC::additional_reference_tangents(
    const int stage,
    const double progress,
    const Eigen::Matrix3d &predictionRotation,
    const Eigen::Matrix4d &predictedParentTransform) const
{
    (void)stage;
    (void)progress;
    (void)predictionRotation;
    (void)predictedParentTransform;
    return Eigen::MatrixXd::Zero(
        ERROR_DIM, stage_control_dimension() - BASE_STAGE_CONTROL_DIMENSION);
}

void
SerialLinkMPCC::configure_additional_stage_inputs(
    const int stage,
    Eigen::VectorXd &lower,
    Eigen::VectorXd &upper,
    Eigen::VectorXd &nominal) const
{
    const int additionalCount =
        stage_control_dimension() - BASE_STAGE_CONTROL_DIMENSION;
    if(additionalCount <= 0)
    {
        return;
    }
    const int offset = stage * stage_control_dimension()
        + BASE_STAGE_CONTROL_DIMENSION;
    lower.segment(offset, additionalCount).setZero();
    upper.segment(offset, additionalCount).setZero();
    nominal.segment(offset, additionalCount).setZero();
}

void
SerialLinkMPCC::reset_additional_virtual_state()
{
}

Eigen::Matrix3d
SerialLinkMPCC::position_error_weight(
    const int stage,
    const Eigen::Vector<double,ERROR_DIM> &pathTangent,
    const Eigen::Matrix3d &defaultWeight,
    const Eigen::Matrix3d &predictionRotation,
    const Eigen::Matrix4d &predictedParentTransform) const
{
    (void)stage;
    (void)pathTangent;
    (void)predictionRotation;
    (void)predictedParentTransform;
    return defaultWeight;
}

void
SerialLinkMPCC::extend_qp_problem(const MpccQpExtensionContext &context,
                                  Eigen::MatrixXd &hessian,
                                  Eigen::VectorXd &gradient,
                                  Eigen::MatrixXd &constraintMatrix,
                                  Eigen::VectorXd &constraintVector,
                                  Eigen::VectorXd &seed)
{
    (void)context;
    (void)hessian;
    (void)gradient;
    (void)constraintMatrix;
    (void)constraintVector;
    (void)seed;
}

void
SerialLinkMPCC::shift_extension_warm_start(
    const MpccQpExtensionContext &context,
    const Eigen::VectorXd &optimum,
    Eigen::VectorXd &shiftedWarmStart)
{
    (void)context;
    (void)optimum;
    (void)shiftedWarmStart;
}

void
SerialLinkMPCC::on_extended_qp_solution(const MpccQpExtensionContext &context,
                                        const Eigen::VectorXd &optimum)
{
    (void)context;
    (void)optimum;
}

Eigen::Vector<double,6>
SerialLinkMPCC::postprocess_base_twist(const Eigen::Vector<double,6> &baseTwist,
                                       const double dt)
{
    (void)dt;
    return baseTwist;
}

void
SerialLinkMPCC::on_twist_resolved(
    const Eigen::Vector<double,6> &commandedBaseTwist,
    const Eigen::Vector<double,6> &realizedBaseTwist,
    const double dt)
{
    (void)commandedBaseTwist;
    (void)realizedBaseTwist;
    (void)dt;
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
