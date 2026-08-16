/**
 * @file    SerialLinkRMPCC.cpp
 * @brief   Riemannian MPCC controller implementation for serial link robot arms.
 */

#include <Control/SerialLinkRMPCC.h>
#include <Control/ProgressSchedule.h>
#include <Control/RmpccPrediction.h>
#include <Control/RmpccProgressConstraints.h>
#include <Math/MathFunctions.h>

#include <Eigen/Geometry>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <vector>

using RobotLibrary::Math::se3_logarithm;
using RobotLibrary::Math::se3_inverse;

namespace RobotLibrary { namespace Control {

namespace {

double clamp_value(const double value, const double lower, const double upper)
{
    return std::max(lower, std::min(upper, value));
}

Eigen::Matrix4d pose_to_matrix(const RobotLibrary::Model::Pose &pose)
{
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    Eigen::Quaterniond q = pose.quaternion();
    q.normalize();
    T.block<3,3>(0,0) = q.toRotationMatrix();
    T.block<3,1>(0,3) = pose.translation();
    return T;
}

RobotLibrary::Model::Pose matrix_to_pose(const Eigen::Matrix4d &T)
{
    Eigen::Quaterniond q(T.block<3,3>(0,0));
    q.normalize();
    return RobotLibrary::Model::Pose(T.block<3,1>(0,3), q);
}

} // anonymous namespace

SerialLinkRMPCC::SerialLinkRMPCC(std::shared_ptr<RobotLibrary::Model::KinematicTree> model,
                                 const std::string &endpointName,
                                 const RobotLibrary::Control::SerialLinkParameters &parameters,
                                 const RmpccParameters &rmpcc)
: SerialLinkVelocityBase(model, endpointName, parameters),
  _rmpcc(rmpcc),
  _deriveProgressRateMin(rmpcc.progressRateMin <= 0.0),
  _deriveProgressRateMax(rmpcc.progressRateMax <= 0.0),
  _qpSolver(parameters.qpsolver)
{
    if(_rmpcc.horizonSteps < 1)
    {
        _rmpcc.horizonSteps = 1;
    }
    if(not std::isfinite(_rmpcc.rtiFiniteDifferenceStep)
       or _rmpcc.rtiFiniteDifferenceStep <= 0.0)
    {
        throw std::invalid_argument(
            "[ERROR] [SERIAL LINK RMPCC] rtiFiniteDifferenceStep must be positive.");
    }
    const auto validateCostScale = [](const double scale, const char *name)
    {
        if(not std::isfinite(scale) or scale < 0.0)
        {
            throw std::invalid_argument(
                std::string("[ERROR] [SERIAL LINK RMPCC] ") + name
                + " must be finite and non-negative.");
        }
    };
    validateCostScale(_rmpcc.runningContourScale, "runningContourScale");
    validateCostScale(_rmpcc.runningLagScale, "runningLagScale");
    validateCostScale(_rmpcc.runningLagTranslationScale, "runningLagTranslationScale");
    validateCostScale(_rmpcc.runningLagRotationScale, "runningLagRotationScale");
    validateCostScale(_rmpcc.terminalLagTranslationScale, "terminalLagTranslationScale");
    validateCostScale(_rmpcc.terminalLagRotationScale, "terminalLagRotationScale");
    validateCostScale(_rmpcc.pathVelocityScale, "pathVelocityScale");
    _lastProgressRate = _rmpcc.progressRateRef;
}

Eigen::VectorXd
SerialLinkRMPCC::track_endpoint_trajectory(const RobotLibrary::Model::Pose &desiredPose,
                                           const Eigen::Vector<double,6>   &desiredVelocity,
                                           const Eigen::Vector<double,6>   &desiredAcceleration)
{
    (void)desiredPose;
    (void)desiredVelocity;
    (void)desiredAcceleration;

    throw std::logic_error(
        "[ERROR] [SERIAL LINK RMPCC] track_endpoint_trajectory(): "
        "Single-pose tracking is not RMPCC. "
        "Call set_trajectory() and then step(dt).");
}

std::string
SerialLinkRMPCC::objective_description() const
{
    std::ostringstream stream;
    stream << "predictor="
           << (_rmpcc.predictorGeometry == RmpccPredictorGeometry::ExactSE3
               ? "exact_se3: E[k+1]=Tref(s[k+1])^-1*Tref(s[k])*E[k]*Exp(dt*u[k])"
               : "additive: e[k+1]=e[k]+dt*(Jr(e[k])^-1*u[k]-g(e[k],s[k])*sdot[k])")
           << "; objective=";
    if(_rmpcc.objectiveGeometry == RmpccObjectiveGeometry::FullScrewSE3)
    {
        stream << "sum ||(I-P_tau)e||_Qc^2+||P_tau*e||_Ql^2"
               << ", P_tau=g*(g^T*M)/(g^T*M*g)";
    }
    else
    {
        stream << "sum ||(I-tp*tp^T)e_p||_Qcp^2+||tp*tp^T*e_p||_Qlp^2"
               << "+||Log(Rref^T*R)||_QR^2; orientation is not projected onto progress";
    }
    stream << "; path_velocity=sum ||u-Ad(E^-1)*tau*sdot||_Rv^2"
           << "; input=sum ||u||_Ru^2; rate=sum ||Delta u||_Rdu^2"
           << "; progress_reward=-q_s*dt*sum(sdot); terminal=last-stage multipliers"
           << "; residual_linearization="
           << (_rmpcc.residualLinearization
                   == RmpccResidualLinearization::FullResidualJacobian
               ? "full_residual_jacobian"
               : "frozen_projector")
           << "; effective_metric_diag=[" << _rmpcc.metric.diagonal().transpose() << "]"
           << "; effective_contour_weight_diag=["
           << _rmpcc.contourWeight.diagonal().transpose() << "]"
           << "; effective_lag_weight_diag=["
           << _rmpcc.lagWeightMatrix.diagonal().transpose() << "]";
    return stream.str();
}

void
SerialLinkRMPCC::set_trajectory(const RobotLibrary::Trajectory::CartesianSpline &trajectory)
{
    const double duration = trajectory.end_time() - trajectory.start_time();
    if(not std::isfinite(duration) or duration <= 0.0)
    {
        throw std::invalid_argument(
            "[ERROR] [SERIAL LINK RMPCC] set_trajectory(): Trajectory duration must be positive.");
    }

    _trajectory = trajectory;
    _trajectorySet = true;
    reset();

    if(_rmpcc.autoProgressRate)
    {
        _rmpcc.progressRateRef = clamp_value(1.0 / duration,
                                             _rmpcc.autoProgressRateMin,
                                             _rmpcc.autoProgressRateMax);
    }
    if(_deriveProgressRateMin)
    {
        _rmpcc.progressRateMin =
            _rmpcc.progressRateRef * _rmpcc.progressRateMinMultiplier;
    }
    if(_deriveProgressRateMax)
    {
        _rmpcc.progressRateMax = _rmpcc.progressRateRef * _rmpcc.progressRateMaxMultiplier;
    }
    _rmpcc.progressRateMax = std::max(_rmpcc.progressRateMax, _rmpcc.progressRateRef);
    _lastProgressRate = _rmpcc.progressRateRef;
}

void
SerialLinkRMPCC::set_trajectory_frame(
    const RobotLibrary::Trajectory::CartesianTrajectoryFrameState &frame)
{
    RobotLibrary::Trajectory::validate_trajectory_frame(frame);
    _trajectoryFrame = frame;
}

void
SerialLinkRMPCC::set_disturbance(const Eigen::Matrix4d &disturbance)
{
    _disturbance = disturbance;
}

void
SerialLinkRMPCC::set_schedule_limit(double scheduleProgressLimit)
{
    _scheduleProgressLimit = clamp_value(scheduleProgressLimit, 0.0, 1.0);
}

void
SerialLinkRMPCC::reset()
{
    _pathProgress = 0.0;
    _lastProgressRate = _rmpcc.progressRateRef;
    _lastBodyTwist.setZero();
    _predictedNextError.setZero();
    _predictionValid = false;
    _warmStart.resize(0);
    _scheduleProgressLimit = 1.0;
    _diagnostics = RmpccDiagnostics();
}

Eigen::Matrix4d
SerialLinkRMPCC::active_frame_transform_in_base() const
{
    return _disturbance * _trajectoryFrame.transformInBase;
}

Eigen::Matrix4d
SerialLinkRMPCC::reference_transform_in_trajectory_frame(double progress)
{
    return pose_to_matrix(_trajectory.pose_at_progress(progress));
}

Eigen::Matrix4d
SerialLinkRMPCC::reference_transform_in_base(double progress)
{
    return active_frame_transform_in_base()
           * reference_transform_in_trajectory_frame(progress);
}

RobotLibrary::Model::Pose
SerialLinkRMPCC::reference_pose(double progress)
{
    return matrix_to_pose(reference_transform_in_base(progress));
}

Eigen::VectorXd
SerialLinkRMPCC::clipped_warm_start(const Eigen::VectorXd &seed,
                                    const Eigen::VectorXd &lower,
                                    const Eigen::VectorXd &upper,
                                    const double dt,
                                    const double remaining,
                                    const double scheduleRemaining) const
{
    Eigen::VectorXd clipped = seed;
    if(clipped.size() != lower.size())
    {
        clipped = Eigen::VectorXd::Zero(lower.size());
    }

    for(int i = 0; i < clipped.size(); ++i)
    {
        clipped(i) = clamp_value(clipped(i), lower(i), upper(i));
    }

    const int N = _rmpcc.horizonSteps;
    if(clipped.size() == 7 * N && N > 0)
    {
        Eigen::VectorXd rates = clipped.segment(6 * N, N);
        rmpcc_clip_progress_rates(
            rates,
            lower.segment(6 * N, N),
            upper.segment(6 * N, N),
            dt,
            remaining + _rmpcc.progressUpperSlack,
            scheduleRemaining);
        clipped.segment(6 * N, N) = rates;
    }

    return clipped;
}

void
SerialLinkRMPCC::solve_rmpcc(const Eigen::Matrix4d &currentTransformInTrajectoryFrame,
                             const double dt)
{
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    _diagnostics.qpStatus = 0.0;
    _diagnostics.fallbackUsed = false;

    const int N = _rmpcc.horizonSteps;
    const int variableDim = 7 * N;
    const int errorDim = TWIST_DIM * N;
    const int progressOffset = TWIST_DIM * N;
    const double remaining = std::max(0.0, 1.0 - _pathProgress);
    const double scheduleRemaining =
        std::max(0.0, _scheduleProgressLimit + _rmpcc.progressScheduleSlack - _pathProgress);

    // Both operands are expressed directly in the active trajectory frame.
    // Their relative transform is invariant to the frame's base transform.
    const Eigen::Matrix4d currentReferenceTransform =
        reference_transform_in_trajectory_frame(_pathProgress);
    const Eigen::Matrix4d relativeTransform =
        se3_inverse(currentReferenceTransform) * currentTransformInTrajectoryFrame;
    const Eigen::Vector<double, TWIST_DIM> e0 = se3_logarithm(relativeTransform);
    _diagnostics.referenceProgress = _pathProgress;
    _diagnostics.modelPredictionResidual =
        _predictionValid ? (e0 - _predictedNextError).norm() : 0.0;
    _diagnostics.se3Error = e0;
    _diagnostics.realizedOneStepErrorNorm = e0.norm();
    const double se3ErrorNorm = e0.norm();
    const double rotationError = e0.tail<3>().norm();
    const Eigen::Vector<double, TWIST_DIM> currentTau =
        _trajectory.tangent_at_progress(_pathProgress, _rmpcc.tangentStep);
    const Eigen::Vector<double, TWIST_DIM> currentErrorTangent =
        rmpcc_error_coordinate_path_tangent(e0, currentTau);
    // Preserve the legacy signed full-screw lag diagnostic even when the
    // running cost uses independent translational/rotational projectors.
    const double currentTauMetric =
        (currentErrorTangent.transpose() * _rmpcc.metric * currentErrorTangent)(0);
    const double currentDenom =
        std::max(currentTauMetric, _rmpcc.hessianRegularization);
    RmpccErrorProjection currentProjection;
    Eigen::Vector<double,TWIST_DIM> currentCostError = e0;
    if(_rmpcc.objectiveGeometry == RmpccObjectiveGeometry::FullScrewSE3)
    {
        currentProjection = rmpcc_error_projection(
            currentErrorTangent, _rmpcc.metric, _rmpcc.runningLagGeometry,
            _rmpcc.hessianRegularization);
    }
    else
    {
        currentCostError = rmpcc_decoupled_error(e0);
        currentProjection = rmpcc_decoupled_error_projection(
            currentTau, _rmpcc.hessianRegularization);
    }
    const Eigen::Matrix<double, TWIST_DIM, TWIST_DIM> &currentLagProjection =
        currentProjection.lag;
    const Eigen::Matrix<double, TWIST_DIM, TWIST_DIM> &currentContourProjection =
        currentProjection.contour;
    const double currentContourError =
        (currentContourProjection * currentCostError).norm();
    const double progressRateMax = _rmpcc.progressRateMax;
    const double progressRateMin = std::min(_rmpcc.progressRateMin, progressRateMax);

    VectorXd lower = VectorXd::Zero(variableDim);
    VectorXd upper = VectorXd::Zero(variableDim);
    const bool relaxForCompletion =
        remaining <= _rmpcc.completionTolerance
        or remaining < static_cast<double>(N) * dt * progressRateMin;
    const VectorXd progressLower = rmpcc_progress_rate_lower_bounds(
        N, progressRateMin, relaxForCompletion, dt, scheduleRemaining);
    const VectorXd fixedProgressRates = reference_schedule_rates(
        _pathProgress, N, dt, 1.0 / _rmpcc.progressRateRef);

    for(int stage = 0; stage < N; ++stage)
    {
        const int uOffset = stage * TWIST_DIM;
        lower.segment<3>(uOffset) = -_rmpcc.linearVelocityMax;
        upper.segment<3>(uOffset) = _rmpcc.linearVelocityMax;
        lower.segment<3>(uOffset + 3) = -_rmpcc.angularVelocityMax;
        upper.segment<3>(uOffset + 3) = _rmpcc.angularVelocityMax;
        lower(progressOffset + stage) = _fixedProgressSchedule
            ? fixedProgressRates(stage) : progressLower(stage);
        upper(progressOffset + stage) = _fixedProgressSchedule
            ? fixedProgressRates(stage) : progressRateMax;
    }

    MatrixXd Bineq = MatrixXd::Zero(2 * variableDim + 2, variableDim);
    VectorXd zineq = VectorXd::Zero(2 * variableDim + 2);
    Bineq.topRows(variableDim).setIdentity();
    zineq.head(variableDim) = upper;
    Bineq.block(variableDim, 0, variableDim, variableDim) = -MatrixXd::Identity(variableDim, variableDim);
    zineq.segment(variableDim, variableDim) = -lower;
    Bineq.block(2 * variableDim, progressOffset, 1, N) =
        rmpcc_completion_progress_row(N, dt).transpose();
    Bineq.block(2 * variableDim + 1, progressOffset, 1, N) =
        rmpcc_schedule_progress_row(N, dt).transpose();
    zineq(2 * variableDim) = remaining + _rmpcc.progressUpperSlack;
    zineq(2 * variableDim + 1) = scheduleRemaining;

    const auto referenceTransform = [this](const double progress)
    {
        return reference_transform_in_trajectory_frame(progress);
    };
    const auto referenceTangent = [this](const double progress)
    {
        return _trajectory.tangent_at_progress(progress, _rmpcc.tangentStep);
    };
    const auto propagate = [&](const RmpccStateVector &state,
                               const RmpccInputVector &input)
    {
        return _rmpcc.predictorGeometry == RmpccPredictorGeometry::ExactSE3
            ? rmpcc_exact_state_step(state, input, dt, referenceTransform)
            : rmpcc_additive_state_step(state, input, dt, referenceTangent);
    };
    const auto linearize = [&](const RmpccStateVector &state,
                               const RmpccInputVector &input)
    {
        return _rmpcc.predictorGeometry == RmpccPredictorGeometry::ExactSE3
            ? rmpcc_linearize_exact_state_step(
                  state, input, dt, _rmpcc.rtiFiniteDifferenceStep,
                  referenceTransform)
            : rmpcc_linearize_additive_state_step(
                  state, input, dt, _rmpcc.rtiFiniteDifferenceStep,
                  referenceTangent);
    };

    if(_warmStart.size() != variableDim)
    {
        _warmStart = VectorXd::Zero(variableDim);
        RmpccStateVector guessState = RmpccStateVector::Zero();
        guessState.head<TWIST_DIM>() = e0;
        guessState(6) = _pathProgress;
        for(int stage = 0; stage < N; ++stage)
        {
            const int uOffset = stage * TWIST_DIM;
            const int sOffset = progressOffset + stage;
            const double rate = clamp_value(_rmpcc.progressRateRef,
                                            lower(sOffset), upper(sOffset));
            const Eigen::Vector<double,TWIST_DIM> tangent =
                _trajectory.tangent_at_progress(guessState(6), _rmpcc.tangentStep);
            Eigen::Vector<double,TWIST_DIM> bodyTwist =
                rmpcc_transport_reference_tangent(
                    guessState.head<TWIST_DIM>(), tangent) * rate;
            for(int component = 0; component < TWIST_DIM; ++component)
            {
                bodyTwist(component) = clamp_value(bodyTwist(component),
                                                   lower(uOffset + component),
                                                   upper(uOffset + component));
            }
            _warmStart.segment<TWIST_DIM>(uOffset) = bodyTwist;
            _warmStart(sOffset) = rate;

            RmpccInputVector guessInput = RmpccInputVector::Zero();
            guessInput.head<TWIST_DIM>() = bodyTwist;
            guessInput(6) = rate;
            guessState = propagate(guessState, guessInput);
        }
    }

    const VectorXd zNominal =
        clipped_warm_start(_warmStart, lower, upper, dt, remaining, scheduleRemaining);

    // Exact nominal rollout followed by one stage-wise RTI/SQP linearisation.
    MatrixXd errorSensitivity = MatrixXd::Zero(errorDim, variableDim);
    VectorXd errorOffset = VectorXd::Zero(errorDim);
    MatrixXd Q = MatrixXd::Zero(errorDim, errorDim);
    MatrixXd fullResidualH = MatrixXd::Zero(variableDim, variableDim);
    VectorXd fullResidualF = VectorXd::Zero(variableDim);
    MatrixXd contourResidualSensitivity = MatrixXd::Zero(errorDim, variableDim);
    VectorXd contourResidualOffset = VectorXd::Zero(errorDim);
    MatrixXd lagResidualSensitivity = MatrixXd::Zero(errorDim, variableDim);
    VectorXd lagResidualOffset = VectorXd::Zero(errorDim);
    Eigen::Matrix<double,7,Eigen::Dynamic> stateSensitivity(7, variableDim);
    stateSensitivity.setZero();

    std::vector<RmpccStateVector> nominalStates(static_cast<size_t>(N + 1));
    std::vector<RmpccInputVector> nominalInputs(static_cast<size_t>(N));
    std::vector<Eigen::Matrix<double,TWIST_DIM,TWIST_DIM>> contourProjections(
        static_cast<size_t>(N));
    std::vector<Eigen::Matrix<double,TWIST_DIM,TWIST_DIM>> lagProjections(
        static_cast<size_t>(N));
    std::vector<Eigen::Matrix<double,TWIST_DIM,TWIST_DIM>> contourCostWeights(
        static_cast<size_t>(N));
    std::vector<Eigen::Matrix<double,TWIST_DIM,TWIST_DIM>> lagCostWeights(
        static_cast<size_t>(N));
    std::vector<Eigen::Vector<double,TWIST_DIM>> transportedTangents(
        static_cast<size_t>(N));
    std::vector<bool> usesCompleteResidualJacobian(static_cast<size_t>(N), false);
    nominalStates.front().setZero();
    nominalStates.front().head<TWIST_DIM>() = e0;
    nominalStates.front()(6) = _pathProgress;

    for(int stage = 0; stage < N; ++stage)
    {
        const int uOffset = stage * TWIST_DIM;
        const int sOffset = progressOffset + stage;
        RmpccInputVector &nominalInput = nominalInputs[static_cast<size_t>(stage)];
        nominalInput.head<TWIST_DIM>() = zNominal.segment<TWIST_DIM>(uOffset);
        nominalInput(6) = zNominal(sOffset);

        const RmpccStageLinearization linearization = linearize(
            nominalStates[static_cast<size_t>(stage)], nominalInput);
        nominalStates[static_cast<size_t>(stage + 1)] = linearization.nominalNext;

        stateSensitivity = linearization.stateJacobian * stateSensitivity;
        for(int component = 0; component < TWIST_DIM; ++component)
        {
            stateSensitivity.col(uOffset + component) +=
                linearization.inputJacobian.col(component);
        }
        stateSensitivity.col(sOffset) += linearization.inputJacobian.col(6);

        const int row = stage * TWIST_DIM;
        const MatrixXd stageSensitivity = stateSensitivity.topRows(TWIST_DIM);

        const double predictedProgress = linearization.nominalNext(6);
        const Eigen::Vector<double,TWIST_DIM> stageReferenceTangent =
            _trajectory.tangent_at_progress(predictedProgress, _rmpcc.tangentStep);
        const Eigen::Vector<double,TWIST_DIM> errorTangent =
            rmpcc_error_coordinate_path_tangent(
                linearization.nominalNext.head<TWIST_DIM>(), stageReferenceTangent);
        // Terminal geometry remains the original full screw projection in all
        // running-lag ablations. This isolates running geometry from terminal
        // semantics just as the leave-one-running-term-out profiles do.
        const RmpccLagGeometry lagGeometry = stage == N - 1
            ? RmpccLagGeometry::FullScrew
            : _rmpcc.runningLagGeometry;
        RmpccErrorProjection stageProjection;
        if(_rmpcc.objectiveGeometry == RmpccObjectiveGeometry::FullScrewSE3)
        {
            stageProjection = rmpcc_error_projection(
                errorTangent, _rmpcc.metric, lagGeometry,
                _rmpcc.hessianRegularization);
            errorSensitivity.block(row, 0, TWIST_DIM, variableDim) = stageSensitivity;
            errorOffset.segment<TWIST_DIM>(row) =
                linearization.nominalNext.head<TWIST_DIM>()
                - stageSensitivity * zNominal;
        }
        else
        {
            const Eigen::Matrix<double,TWIST_DIM,TWIST_DIM> errorMap =
                rmpcc_decoupled_error_jacobian(
                    linearization.nominalNext.head<TWIST_DIM>(),
                    _rmpcc.rtiFiniteDifferenceStep);
            const MatrixXd mappedSensitivity = errorMap * stageSensitivity;
            errorSensitivity.block(row, 0, TWIST_DIM, variableDim) = mappedSensitivity;
            errorOffset.segment<TWIST_DIM>(row) =
                rmpcc_decoupled_error(linearization.nominalNext.head<TWIST_DIM>())
                - mappedSensitivity * zNominal;
            stageProjection = rmpcc_decoupled_error_projection(
                stageReferenceTangent, _rmpcc.hessianRegularization);
        }
        const Eigen::Matrix<double,TWIST_DIM,TWIST_DIM> &lagProjection =
            stageProjection.lag;
        const Eigen::Matrix<double,TWIST_DIM,TWIST_DIM> &contourProjection =
            stageProjection.contour;
        Eigen::Matrix<double,TWIST_DIM,TWIST_DIM> stageContourWeight =
            _rmpcc.contourWeight;
        Eigen::Matrix<double,TWIST_DIM,TWIST_DIM> stageLagWeight =
            _rmpcc.lagWeightMatrix;
        if(stage == N - 1)
        {
            const double positionScale =
                std::sqrt(_rmpcc.terminalPositionMultiplier);
            const double rotationScale =
                std::sqrt(_rmpcc.terminalRotationMultiplier);
            Eigen::Matrix<double,TWIST_DIM,TWIST_DIM> terminalScale =
                Eigen::Matrix<double,TWIST_DIM,TWIST_DIM>::Identity();
            terminalScale.diagonal().head<3>().setConstant(positionScale);
            terminalScale.diagonal().tail<3>().setConstant(rotationScale);
            stageContourWeight = terminalScale * stageContourWeight * terminalScale;
            stageLagWeight *= _rmpcc.terminalLagMultiplier;
            stageLagWeight = rmpcc_component_scaled_weight(
                stageLagWeight,
                _rmpcc.terminalLagTranslationScale,
                _rmpcc.terminalLagRotationScale);
        }
        else
        {
            // Running-cost ablations deliberately leave the terminal objective
            // unchanged so each experiment removes exactly one running term.
            stageContourWeight *= _rmpcc.runningContourScale;
            stageLagWeight *= _rmpcc.runningLagScale;
            stageLagWeight = rmpcc_component_scaled_weight(
                stageLagWeight,
                _rmpcc.runningLagTranslationScale,
                _rmpcc.runningLagRotationScale);
        }
        Eigen::Matrix<double,TWIST_DIM,TWIST_DIM> stageWeight;
        if(_rmpcc.objectiveGeometry == RmpccObjectiveGeometry::FullScrewSE3)
        {
            stageWeight =
                contourProjection.transpose() * stageContourWeight * contourProjection
                + lagProjection.transpose() * stageLagWeight * lagProjection;
        }
        else
        {
            stageWeight = rmpcc_decoupled_cost_weight(
                stageReferenceTangent, stageContourWeight, stageLagWeight,
                _rmpcc.hessianRegularization);
        }
        if(stage == N - 1)
        {
            // Backward-compatible global multiplier; new configurations should
            // normally leave it at one and tune the three terminal terms above.
            stageWeight *= _rmpcc.terminalMultiplier;
            stageContourWeight *= _rmpcc.terminalMultiplier;
            stageLagWeight *= _rmpcc.terminalMultiplier;
        }
        contourProjections[static_cast<size_t>(stage)] = contourProjection;
        lagProjections[static_cast<size_t>(stage)] = lagProjection;
        contourCostWeights[static_cast<size_t>(stage)] = stageContourWeight;
        lagCostWeights[static_cast<size_t>(stage)] = stageLagWeight;
        const bool useFullResidualJacobian =
            _rmpcc.residualLinearization
                == RmpccResidualLinearization::FullResidualJacobian
            and (_rmpcc.objectiveGeometry
                    == RmpccObjectiveGeometry::DecoupledCartesianSO3
                 or lagGeometry == RmpccLagGeometry::FullScrew);
        if(useFullResidualJacobian)
        {
            const RmpccFullScrewResidualLinearization residualLinearization =
                _rmpcc.objectiveGeometry == RmpccObjectiveGeometry::FullScrewSE3
                ? rmpcc_linearize_full_screw_residuals(
                      linearization.nominalNext,
                      _rmpcc.metric,
                      _rmpcc.hessianRegularization,
                      _rmpcc.rtiFiniteDifferenceStep,
                      referenceTangent)
                : rmpcc_linearize_decoupled_residuals(
                      linearization.nominalNext,
                      _rmpcc.hessianRegularization,
                      _rmpcc.rtiFiniteDifferenceStep,
                      referenceTangent);
            const MatrixXd contourJacobian =
                residualLinearization.contourJacobian * stateSensitivity;
            const MatrixXd lagJacobian =
                residualLinearization.lagJacobian * stateSensitivity;
            const Eigen::Vector<double,TWIST_DIM> contourOffset =
                residualLinearization.residual.contour
                - contourJacobian * zNominal;
            const Eigen::Vector<double,TWIST_DIM> lagOffset =
                residualLinearization.residual.lag
                - lagJacobian * zNominal;

            fullResidualH +=
                2.0 * contourJacobian.transpose()
                    * stageContourWeight * contourJacobian
                + 2.0 * lagJacobian.transpose()
                    * stageLagWeight * lagJacobian;
            fullResidualF +=
                2.0 * contourJacobian.transpose()
                    * stageContourWeight * contourOffset
                + 2.0 * lagJacobian.transpose()
                    * stageLagWeight * lagOffset;
            contourResidualSensitivity.block(row, 0, TWIST_DIM, variableDim) =
                contourJacobian;
            contourResidualOffset.segment<TWIST_DIM>(row) = contourOffset;
            lagResidualSensitivity.block(row, 0, TWIST_DIM, variableDim) =
                lagJacobian;
            lagResidualOffset.segment<TWIST_DIM>(row) = lagOffset;
            usesCompleteResidualJacobian[static_cast<size_t>(stage)] = true;
        }
        else
        {
            Q.block<TWIST_DIM,TWIST_DIM>(row, row) = stageWeight;
            if(_rmpcc.objectiveGeometry == RmpccObjectiveGeometry::FullScrewSE3)
            {
                contourResidualSensitivity.block(row, 0, TWIST_DIM, variableDim) =
                    contourProjection * stageSensitivity;
                contourResidualOffset.segment<TWIST_DIM>(row) =
                    contourProjection * errorOffset.segment<TWIST_DIM>(row);
                lagResidualSensitivity.block(row, 0, TWIST_DIM, variableDim) =
                    lagProjection * stageSensitivity;
                lagResidualOffset.segment<TWIST_DIM>(row) =
                    lagProjection * errorOffset.segment<TWIST_DIM>(row);
            }
        }
    }

    MatrixXd H = 2.0 * errorSensitivity.transpose() * Q * errorSensitivity
                 + fullResidualH;
    VectorXd f = 2.0 * errorSensitivity.transpose() * Q * errorOffset
                 + fullResidualF;

    for(int stage = 0; stage < N; ++stage)
    {
        const int uOffset = stage * TWIST_DIM;
        const int sOffset = progressOffset + stage;
        const RmpccStateVector &stageState = nominalStates[static_cast<size_t>(stage)];
        const Eigen::Vector<double,TWIST_DIM> tangent =
            _trajectory.tangent_at_progress(stageState(6), _rmpcc.tangentStep);
        const Eigen::Vector<double,TWIST_DIM> transportedTangent =
            rmpcc_transport_reference_tangent(
                stageState.head<TWIST_DIM>(), tangent);
        transportedTangents[static_cast<size_t>(stage)] = transportedTangent;
        const Eigen::Matrix<double,TWIST_DIM,TWIST_DIM> pathVelocityWeight =
            _rmpcc.pathVelocityScale * _rmpcc.pathVelocityWeight;

        H.block<TWIST_DIM,TWIST_DIM>(uOffset, uOffset) +=
            2.0 * (_rmpcc.controlWeight + pathVelocityWeight);
        H.block<TWIST_DIM,1>(uOffset, sOffset) +=
            -2.0 * pathVelocityWeight * transportedTangent;
        H.block<1,TWIST_DIM>(sOffset, uOffset) +=
            -2.0 * transportedTangent.transpose() * pathVelocityWeight;
        H(sOffset, sOffset) +=
            2.0 * (transportedTangent.transpose()
                   * pathVelocityWeight * transportedTangent)(0);

        if(_rmpcc.progressRateWeight > 0.0)
        {
            H(sOffset, sOffset) += 2.0 * _rmpcc.progressRateWeight;
            f(sOffset) += -2.0 * _rmpcc.progressRateWeight * _rmpcc.progressRateRef;
        }
        f(sOffset) += -_rmpcc.progressReward * dt;

        if(stage == 0)
        {
            H.block<TWIST_DIM,TWIST_DIM>(uOffset, uOffset) +=
                2.0 * _rmpcc.controlRateWeight;
            f.segment<TWIST_DIM>(uOffset) +=
                -2.0 * _rmpcc.controlRateWeight * _lastBodyTwist;
            H(sOffset, sOffset) += 2.0 * _rmpcc.progressRateSmoothWeight;
            f(sOffset) += -2.0 * _rmpcc.progressRateSmoothWeight * _lastProgressRate;
        }
        else
        {
            const int previousUOffset = (stage - 1) * TWIST_DIM;
            const int previousSOffset = progressOffset + stage - 1;
            H.block<TWIST_DIM,TWIST_DIM>(uOffset, uOffset) +=
                2.0 * _rmpcc.controlRateWeight;
            H.block<TWIST_DIM,TWIST_DIM>(previousUOffset, previousUOffset) +=
                2.0 * _rmpcc.controlRateWeight;
            H.block<TWIST_DIM,TWIST_DIM>(uOffset, previousUOffset) +=
                -2.0 * _rmpcc.controlRateWeight;
            H.block<TWIST_DIM,TWIST_DIM>(previousUOffset, uOffset) +=
                -2.0 * _rmpcc.controlRateWeight;
            H(sOffset, sOffset) += 2.0 * _rmpcc.progressRateSmoothWeight;
            H(previousSOffset, previousSOffset) +=
                2.0 * _rmpcc.progressRateSmoothWeight;
            H(sOffset, previousSOffset) += -2.0 * _rmpcc.progressRateSmoothWeight;
            H(previousSOffset, sOffset) += -2.0 * _rmpcc.progressRateSmoothWeight;
        }
    }

    H += _rmpcc.hessianRegularization * MatrixXd::Identity(variableDim, variableDim);
    H = 0.5 * (H + H.transpose());

    const auto solveStart = std::chrono::steady_clock::now();
    VectorXd zOpt = _qpSolver.solve(H, f, Bineq, zineq, zNominal);
    _diagnostics.qpSolveTimeSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - solveStart).count();
    if(zOpt.size() != variableDim or not zOpt.allFinite())
    {
        throw std::runtime_error("[ERROR] [SERIAL LINK RMPCC] solve_rmpcc(): QP returned an invalid solution.");
    }
    const double constraintViolation = (Bineq * zOpt - zineq).maxCoeff();
    if(not std::isfinite(constraintViolation) or constraintViolation > 1e-6)
    {
        throw std::runtime_error(
            "[ERROR] [SERIAL LINK RMPCC] solve_rmpcc(): QP returned an infeasible solution (maximum violation "
            + std::to_string(constraintViolation) + ").");
    }
    const VectorXd qpGradient = H * zOpt + f;
    _diagnostics.qpFirstTwistGradientNorm = qpGradient.head<TWIST_DIM>().norm();
    _diagnostics.qpStepNorm = (zOpt - zNominal).norm();
    _diagnostics.activeConstraintCount = 0.0;
    const Eigen::VectorXd constraintSlack = zineq - Bineq * zOpt;
    for(int i = 0; i < constraintSlack.size(); ++i)
    {
        _diagnostics.activeConstraintCount += constraintSlack(i) <= 1e-6 ? 1.0 : 0.0;
    }

    zOpt = clipped_warm_start(zOpt, lower, upper, dt, remaining, scheduleRemaining);

    _warmStart = zOpt;
    if(_warmStart.size() == variableDim)
    {
        for(int stage = 0; stage < N - 1; ++stage)
        {
            _warmStart.segment<TWIST_DIM>(stage * TWIST_DIM) =
                zOpt.segment<TWIST_DIM>((stage + 1) * TWIST_DIM);
            _warmStart(progressOffset + stage) = zOpt(progressOffset + stage + 1);
        }
        _warmStart.segment<TWIST_DIM>((N - 1) * TWIST_DIM) =
            zOpt.segment<TWIST_DIM>((N - 1) * TWIST_DIM);
        _warmStart(progressOffset + N - 1) = zOpt(progressOffset + N - 1);
    }

    _diagnostics.bodyTwist = zOpt.head<TWIST_DIM>();
    _diagnostics.progressRate = zOpt(progressOffset);
    _diagnostics.se3ErrorNorm = se3ErrorNorm;
    _diagnostics.translationError = e0.head<3>().norm();
    _diagnostics.rotationError = rotationError;
    _diagnostics.qpStatus = 1.0;
    _diagnostics.referenceLinearSpeed = (currentTau * _rmpcc.progressRateRef).head<3>().norm();
    _diagnostics.referenceAngularSpeed = (currentTau * _rmpcc.progressRateRef).tail<3>().norm();
    _diagnostics.feedforwardBodyTwist =
        rmpcc_transport_reference_tangent(e0, currentTau)
        * _diagnostics.progressRate;
    _diagnostics.correctionBodyTwist =
        _diagnostics.bodyTwist - _diagnostics.feedforwardBodyTwist;
    _diagnostics.linearVelocityLimitActive = 0.0;
    _diagnostics.angularVelocityLimitActive = 0.0;
    for(int component = 0; component < 3; ++component)
    {
        _diagnostics.linearVelocityLimitActive = std::max(
            _diagnostics.linearVelocityLimitActive,
            std::abs(_diagnostics.bodyTwist(component))
                >= _rmpcc.linearVelocityMax(component) - 1e-6 ? 1.0 : 0.0);
        _diagnostics.angularVelocityLimitActive = std::max(
            _diagnostics.angularVelocityLimitActive,
            std::abs(_diagnostics.bodyTwist(component + 3))
                >= _rmpcc.angularVelocityMax(component) - 1e-6 ? 1.0 : 0.0);
    }
    _diagnostics.maxAngularComponent =
        _diagnostics.bodyTwist.tail<3>().cwiseAbs().maxCoeff();
    const Eigen::Vector<double,TWIST_DIM> pathVelocityResidual =
        _diagnostics.bodyTwist
        - rmpcc_transport_reference_tangent(e0, currentTau)
          * _diagnostics.progressRate;
    _diagnostics.pathVelocityLinearResidual =
        pathVelocityResidual.head<3>().norm();
    _diagnostics.pathVelocityAngularResidual =
        pathVelocityResidual.tail<3>().norm();
    const Eigen::Vector<double,TWIST_DIM> currentLagError =
        currentLagProjection * currentCostError;
    const Eigen::Matrix<double, 1, TWIST_DIM> lagRow =
        (currentErrorTangent.transpose() * _rmpcc.metric) / std::sqrt(currentDenom);
    _diagnostics.lagError = (lagRow * e0)(0);
    _diagnostics.lagErrorNorm = currentLagError.norm();
    _diagnostics.lagTranslationErrorNorm = currentLagError.head<3>().norm();
    _diagnostics.lagRotationErrorNorm = currentLagError.tail<3>().norm();
    _diagnostics.contourError = currentContourError;

    _diagnostics.runningContourCost = 0.0;
    _diagnostics.runningLagCost = 0.0;
    _diagnostics.runningPathVelocityCost = 0.0;
    _diagnostics.terminalContourCost = 0.0;
    _diagnostics.terminalLagCost = 0.0;
    _diagnostics.terminalLagTranslationCost = 0.0;
    _diagnostics.terminalLagRotationCost = 0.0;
    _diagnostics.runningLagTranslationCost = 0.0;
    _diagnostics.runningLagRotationCost = 0.0;
    const VectorXd predictedErrors = errorSensitivity * zOpt + errorOffset;
    const VectorXd predictedContourResiduals =
        contourResidualSensitivity * zOpt + contourResidualOffset;
    const VectorXd predictedLagResiduals =
        lagResidualSensitivity * zOpt + lagResidualOffset;
    const Eigen::Matrix<double,TWIST_DIM,TWIST_DIM> pathVelocityWeight =
        _rmpcc.pathVelocityScale * _rmpcc.pathVelocityWeight;
    for(int stage = 0; stage < N; ++stage)
    {
        const int row = stage * TWIST_DIM;
        const int uOffset = stage * TWIST_DIM;
        const int sOffset = progressOffset + stage;
        const Eigen::Vector<double,TWIST_DIM> predictedError =
            predictedErrors.segment<TWIST_DIM>(row);
        Eigen::Vector<double,TWIST_DIM> contourError;
        Eigen::Vector<double,TWIST_DIM> lagError;
        if(_rmpcc.objectiveGeometry == RmpccObjectiveGeometry::FullScrewSE3
           or usesCompleteResidualJacobian[static_cast<size_t>(stage)])
        {
            contourError = predictedContourResiduals.segment<TWIST_DIM>(row);
            lagError = predictedLagResiduals.segment<TWIST_DIM>(row);
        }
        else
        {
            contourError =
                contourProjections[static_cast<size_t>(stage)] * predictedError;
            lagError = lagProjections[static_cast<size_t>(stage)] * predictedError;
        }
        const double contourCost =
            (contourError.transpose()
             * contourCostWeights[static_cast<size_t>(stage)] * contourError)(0);
        const double lagCost =
            (lagError.transpose()
             * lagCostWeights[static_cast<size_t>(stage)] * lagError)(0);
        if(stage == N - 1)
        {
            _diagnostics.terminalContourCost = contourCost;
            _diagnostics.terminalLagCost = lagCost;
            _diagnostics.terminalLagTranslationCost =
                (lagError.head<3>().transpose()
                 * lagCostWeights[static_cast<size_t>(stage)].block<3,3>(0,0)
                 * lagError.head<3>())(0);
            _diagnostics.terminalLagRotationCost =
                (lagError.tail<3>().transpose()
                 * lagCostWeights[static_cast<size_t>(stage)].block<3,3>(3,3)
                 * lagError.tail<3>())(0);
        }
        else
        {
            _diagnostics.runningContourCost += contourCost;
            _diagnostics.runningLagCost += lagCost;
            _diagnostics.runningLagTranslationCost +=
                (lagError.head<3>().transpose()
                 * lagCostWeights[static_cast<size_t>(stage)].block<3,3>(0,0)
                 * lagError.head<3>())(0);
            _diagnostics.runningLagRotationCost +=
                (lagError.tail<3>().transpose()
                 * lagCostWeights[static_cast<size_t>(stage)].block<3,3>(3,3)
                 * lagError.tail<3>())(0);
        }

        const Eigen::Vector<double,TWIST_DIM> residual =
            zOpt.segment<TWIST_DIM>(uOffset)
            - transportedTangents[static_cast<size_t>(stage)] * zOpt(sOffset);
        _diagnostics.runningPathVelocityCost +=
            (residual.transpose() * pathVelocityWeight * residual)(0);
    }
}

Eigen::VectorXd
SerialLinkRMPCC::step(const double dt)
{
    if(not _trajectorySet)
    {
        throw std::runtime_error("[ERROR] [SERIAL LINK RMPCC] step(): No trajectory set. Call set_trajectory() first.");
    }
    if(not std::isfinite(dt) or dt <= 0.0)
    {
        throw std::invalid_argument("[ERROR] [SERIAL LINK RMPCC] step(): dt must be positive.");
    }
    update();

    const RobotLibrary::Model::Pose currentPose = endpoint_pose();
    const Eigen::Matrix4d currentTransformInBase = pose_to_matrix(currentPose);
    const Eigen::Matrix4d currentTransformInTrajectoryFrame =
        se3_inverse(active_frame_transform_in_base()) * currentTransformInBase;

    solve_rmpcc(currentTransformInTrajectoryFrame, dt);
    _diagnostics.effectiveLoopFrequency = 1.0 / dt;

    // The model Jacobian maps qdot to endpoint point velocity [p_dot; omega],
    // not to the screw-theory spatial twist [v; omega]. Rotate both components
    // into the base frame without the SE(3) adjoint p x omega term.
    const Eigen::Matrix3d currentRotation = currentTransformInBase.block<3,3>(0,0);
    Eigen::Vector<double, TWIST_DIM> baseTwist;
    baseTwist.head<3>() = currentRotation * _diagnostics.bodyTwist.head<3>();
    baseTwist.tail<3>() = currentRotation * _diagnostics.bodyTwist.tail<3>();

    const Eigen::VectorXd jointCommand = resolve_endpoint_twist(baseTwist);
    _diagnostics.jointVelocityLimitActive = 0.0;
    for(int joint = 0; joint < jointCommand.size(); ++joint)
    {
        const RobotLibrary::Model::Limits limits =
            compute_control_limits(static_cast<unsigned int>(joint));
        if(jointCommand(joint) <= limits.lower + 1e-3
           or jointCommand(joint) >= limits.upper - 1e-3)
        {
            _diagnostics.jointVelocityLimitActive = 1.0;
            break;
        }
    }
    const Eigen::Vector<double, TWIST_DIM> realizedBaseTwist = _jacobianMatrix * jointCommand;
    Eigen::Vector<double, TWIST_DIM> realizedBodyTwist;
    realizedBodyTwist.head<3>() = currentRotation.transpose() * realizedBaseTwist.head<3>();
    realizedBodyTwist.tail<3>() = currentRotation.transpose() * realizedBaseTwist.tail<3>();
    _diagnostics.twistRealizationError =
        (realizedBodyTwist - _diagnostics.bodyTwist).norm();

    RmpccStateVector realizedState = RmpccStateVector::Zero();
    realizedState.head<TWIST_DIM>() = _diagnostics.se3Error;
    realizedState(6) = _pathProgress;
    RmpccInputVector realizedInput = RmpccInputVector::Zero();
    realizedInput.head<TWIST_DIM>() = realizedBodyTwist;
    realizedInput(6) = _diagnostics.progressRate;
    const auto referenceTransform = [this](const double progress)
    {
        return reference_transform_in_trajectory_frame(progress);
    };
    if(_rmpcc.predictorGeometry == RmpccPredictorGeometry::ExactSE3)
    {
        _predictedNextError =
            rmpcc_exact_state_step(realizedState, realizedInput, dt,
                                   referenceTransform).head<TWIST_DIM>();
    }
    else
    {
        const auto referenceTangent = [this](const double progress)
        {
            return _trajectory.tangent_at_progress(progress, _rmpcc.tangentStep);
        };
        _predictedNextError =
            rmpcc_additive_state_step(realizedState, realizedInput, dt,
                                      referenceTangent).head<TWIST_DIM>();
    }
    _diagnostics.predictedNextErrorNorm = _predictedNextError.norm();
    _predictionValid = true;

    _lastBodyTwist = _diagnostics.bodyTwist;
    _lastProgressRate = _diagnostics.progressRate;
    _pathProgress = clamp_value(_pathProgress + dt * _diagnostics.progressRate, 0.0, 1.0);
    _diagnostics.pathProgress = _pathProgress;

    return jointCommand;
}

Eigen::VectorXd
SerialLinkRMPCC::step(const double dt, const double estimatedProgress)
{
    if(not std::isfinite(estimatedProgress))
        throw std::invalid_argument("[ERROR] [SERIAL LINK RMPCC] step(): estimatedProgress must be finite.");
    _pathProgress = clamp_value(estimatedProgress, 0.0, 1.0);
    return step(dt);
}

} } // namespace
