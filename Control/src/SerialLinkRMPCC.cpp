/**
 * @file    SerialLinkRMPCC.cpp
 * @brief   Riemannian MPCC controller implementation for serial link robot arms.
 */

#include <Control/SerialLinkRMPCC.h>
#include <Math/MathFunctions.h>

#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

using RobotLibrary::Math::se3_logarithm;
using RobotLibrary::Math::se3_inverse;
using RobotLibrary::Math::so3_logarithm;

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

// World-frame pose error using axis-angle (so3 logarithm) for orientation, consistent
// with the se3 logarithm convention used inside the QP.
Eigen::Vector<double,6> pose_feedback_error(const RobotLibrary::Model::Pose &current,
                                            const RobotLibrary::Model::Pose &desired)
{
    Eigen::Vector<double,6> error = Eigen::Vector<double,6>::Zero();
    error.head<3>() = desired.translation() - current.translation();
    const Eigen::Matrix3d R_err =
        desired.quaternion().toRotationMatrix() * current.quaternion().inverse().toRotationMatrix();
    error.tail<3>() = so3_logarithm(R_err);
    return error;
}

} // anonymous namespace

SerialLinkRMPCC::SerialLinkRMPCC(std::shared_ptr<RobotLibrary::Model::KinematicTree> model,
                                 const std::string &endpointName,
                                 const RobotLibrary::Control::SerialLinkParameters &parameters,
                                 const RmpccParameters &rmpcc)
: SerialLinkVelocityBase(model, endpointName, parameters),
  _rmpcc(rmpcc),
  _qpSolver(parameters.qpsolver)
{
    if(_rmpcc.horizonSteps < 1)
    {
        _rmpcc.horizonSteps = 1;
    }
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
        "Call set_trajectory() and then step(dt, estimatedProgress).");
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
    if(_rmpcc.progressRateMax <= 0.0)
    {
        _rmpcc.progressRateMax = _rmpcc.progressRateRef * _rmpcc.progressRateMaxMultiplier;
    }
    _rmpcc.progressRateMax = std::max(_rmpcc.progressRateMax, _rmpcc.progressRateRef);
    _lastProgressRate = _rmpcc.progressRateRef;
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
    _warmStart.resize(0);
    _scheduleProgressLimit = 1.0;
    _diagnostics = RmpccDiagnostics();
}

Eigen::Matrix4d
SerialLinkRMPCC::reference_transform(double progress)
{
    // Nominal path geometry comes from the trajectory; the rigid disturbance is a
    // controller-level concern applied on top (it cancels in the tangent, only shifts e0).
    return _disturbance * pose_to_matrix(_trajectory.pose_at_progress(progress));
}

RobotLibrary::Model::Pose
SerialLinkRMPCC::reference_pose(double progress)
{
    return matrix_to_pose(reference_transform(progress));
}

Eigen::Vector<double, SerialLinkRMPCC::TWIST_DIM>
SerialLinkRMPCC::body_twist_reference_at_progress(double progress)
{
    return _trajectory.tangent_at_progress(progress, _rmpcc.tangentStep) * _rmpcc.progressRateRef;
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
        const double maxTotalProgress =
            std::min(remaining + _rmpcc.progressUpperSlack, scheduleRemaining);
        double totalProgress = 0.0;
        for(int stage = 0; stage < N; ++stage)
        {
            totalProgress += dt * clipped(6 * N + stage);
        }

        if(totalProgress > maxTotalProgress)
        {
            const double feasibleRate =
                std::max(0.0, maxTotalProgress) / (static_cast<double>(N) * std::max(dt, 1e-9));
            for(int stage = 0; stage < N; ++stage)
            {
                clipped(6 * N + stage) =
                    clamp_value(feasibleRate, lower(6 * N + stage), upper(6 * N + stage));
            }
        }
    }

    return clipped;
}

void
SerialLinkRMPCC::solve_rmpcc(const Eigen::Matrix4d &currentTransform, const double dt)
{
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    const int N = _rmpcc.horizonSteps;
    const int variableDim = 7 * N;
    const int errorDim = TWIST_DIM * N;
    const int progressOffset = TWIST_DIM * N;
    const double remaining = std::max(0.0, 1.0 - _pathProgress);
    const double scheduleRemaining =
        std::max(0.0, _scheduleProgressLimit + _rmpcc.progressScheduleSlack - _pathProgress);

    const Eigen::Matrix4d referenceTransform = reference_transform(_pathProgress);
    const Eigen::Vector<double, TWIST_DIM> e0 = se3_logarithm(se3_inverse(referenceTransform) * currentTransform);
    const double se3ErrorNorm = e0.norm();
    const double rotationError = e0.tail<3>().norm();
    const Eigen::Vector<double, TWIST_DIM> currentTau = _trajectory.tangent_at_progress(_pathProgress, _rmpcc.tangentStep);
    const double currentTauMetric = (currentTau.transpose() * _rmpcc.metric * currentTau)(0);
    const double currentDenom = std::max(currentTauMetric, _rmpcc.hessianRegularization);
    const Eigen::Matrix<double, TWIST_DIM, TWIST_DIM> currentLagProjection =
        currentTau * (currentTau.transpose() * _rmpcc.metric) / currentDenom;
    const Eigen::Matrix<double, TWIST_DIM, TWIST_DIM> currentContourProjection =
        Eigen::Matrix<double, TWIST_DIM, TWIST_DIM>::Identity() - currentLagProjection;
    const double currentContourError = (currentContourProjection * e0).norm();
    const double progressRateMax = _rmpcc.progressRateMax;
    const double progressRateMin = std::min(_rmpcc.progressRateMin, progressRateMax);

    MatrixXd A = MatrixXd::Zero(errorDim, variableDim);
    VectorXd b = VectorXd::Zero(errorDim);
    MatrixXd Q = MatrixXd::Zero(errorDim, errorDim);

    std::vector<double> predictedProgress(static_cast<size_t>(N), _pathProgress);
    double accumulatedProgress = _pathProgress;
    for(int stage = 0; stage < N; ++stage)
    {
        double sdotGuess = _rmpcc.progressRateRef;
        if(_warmStart.size() == variableDim)
        {
            sdotGuess = _warmStart(progressOffset + stage);
        }
        accumulatedProgress += dt * clamp_value(sdotGuess, progressRateMin, progressRateMax);
        predictedProgress[static_cast<size_t>(stage)] = clamp_value(accumulatedProgress, 0.0, 1.0);
    }

    for(int stage = 0; stage < N; ++stage)
    {
        const int row = stage * TWIST_DIM;
        b.segment<TWIST_DIM>(row) = e0;

        const double s = predictedProgress[static_cast<size_t>(stage)];
        const Eigen::Vector<double, TWIST_DIM> tau = _trajectory.tangent_at_progress(s, _rmpcc.tangentStep);
        const double tauMetric = (tau.transpose() * _rmpcc.metric * tau)(0);
        const double denom = std::max(tauMetric, _rmpcc.hessianRegularization);
        const Eigen::Matrix<double, 1, TWIST_DIM> lagRow =
            (tau.transpose() * _rmpcc.metric) / std::sqrt(denom);
        const Eigen::Matrix<double, TWIST_DIM, TWIST_DIM> lagProjection =
            tau * (tau.transpose() * _rmpcc.metric) / denom;
        const Eigen::Matrix<double, TWIST_DIM, TWIST_DIM> contourProjection =
            Eigen::Matrix<double, TWIST_DIM, TWIST_DIM>::Identity() - lagProjection;
        Eigen::Matrix<double, TWIST_DIM, TWIST_DIM> stageWeight =
            contourProjection.transpose() * _rmpcc.contourWeight * contourProjection
            + _rmpcc.lagWeight * lagRow.transpose() * lagRow;
        if(stage == N - 1)
        {
            stageWeight *= _rmpcc.terminalMultiplier;
        }
        Q.block<TWIST_DIM, TWIST_DIM>(row, row) = stageWeight;

        for(int input = 0; input <= stage; ++input)
        {
            const double inputProgress = predictedProgress[static_cast<size_t>(input)];
            const Eigen::Vector<double, TWIST_DIM> inputTau = _trajectory.tangent_at_progress(inputProgress, _rmpcc.tangentStep);
            A.block<TWIST_DIM, TWIST_DIM>(row, input * TWIST_DIM) +=
                dt * Eigen::Matrix<double, TWIST_DIM, TWIST_DIM>::Identity();
            A.block<TWIST_DIM, 1>(row, progressOffset + input) += -dt * inputTau;
        }
    }

    MatrixXd H = 2.0 * A.transpose() * Q * A;
    VectorXd f = 2.0 * A.transpose() * Q * b;

    for(int stage = 0; stage < N; ++stage)
    {
        const int uOffset = stage * TWIST_DIM;
        const int sOffset = progressOffset + stage;
        const double feedforwardProgress = predictedProgress[static_cast<size_t>(stage)];
        const Eigen::Vector<double, TWIST_DIM> tau =
            _trajectory.tangent_at_progress(feedforwardProgress, _rmpcc.tangentStep);

        H.block<TWIST_DIM, TWIST_DIM>(uOffset, uOffset) += 2.0 * _rmpcc.controlWeight;
        H.block<TWIST_DIM, TWIST_DIM>(uOffset, uOffset) += 2.0 * _rmpcc.pathVelocityWeight;
        H.block<TWIST_DIM, 1>(uOffset, sOffset) +=
            -2.0 * _rmpcc.pathVelocityWeight * tau;
        H.block<1, TWIST_DIM>(sOffset, uOffset) +=
            -2.0 * tau.transpose() * _rmpcc.pathVelocityWeight;
        H(sOffset, sOffset) +=
            2.0 * (tau.transpose() * _rmpcc.pathVelocityWeight * tau)(0);

        if(_rmpcc.progressRateWeight > 0.0)
        {
            H(sOffset, sOffset) += 2.0 * _rmpcc.progressRateWeight;
            f(sOffset) += -2.0 * _rmpcc.progressRateWeight * _rmpcc.progressRateRef;
        }
        f(sOffset) += -_rmpcc.progressReward;

        if(stage == 0)
        {
            H.block<TWIST_DIM, TWIST_DIM>(uOffset, uOffset) += 2.0 * _rmpcc.controlRateWeight;
            f.segment<TWIST_DIM>(uOffset) += -2.0 * _rmpcc.controlRateWeight * _lastBodyTwist;
        }
        else
        {
            const int previousOffset = (stage - 1) * TWIST_DIM;
            H.block<TWIST_DIM, TWIST_DIM>(uOffset, uOffset) += 2.0 * _rmpcc.controlRateWeight;
            H.block<TWIST_DIM, TWIST_DIM>(previousOffset, previousOffset) += 2.0 * _rmpcc.controlRateWeight;
            H.block<TWIST_DIM, TWIST_DIM>(uOffset, previousOffset) += -2.0 * _rmpcc.controlRateWeight;
            H.block<TWIST_DIM, TWIST_DIM>(previousOffset, uOffset) += -2.0 * _rmpcc.controlRateWeight;
        }

        if(stage == 0)
        {
            H(sOffset, sOffset) += 2.0 * _rmpcc.progressRateSmoothWeight;
            f(sOffset) += -2.0 * _rmpcc.progressRateSmoothWeight * _lastProgressRate;
        }
        else
        {
            const int previousSOffset = progressOffset + stage - 1;
            H(sOffset, sOffset) += 2.0 * _rmpcc.progressRateSmoothWeight;
            H(previousSOffset, previousSOffset) += 2.0 * _rmpcc.progressRateSmoothWeight;
            H(sOffset, previousSOffset) += -2.0 * _rmpcc.progressRateSmoothWeight;
            H(previousSOffset, sOffset) += -2.0 * _rmpcc.progressRateSmoothWeight;
        }
    }

    H += _rmpcc.hessianRegularization * MatrixXd::Identity(variableDim, variableDim);
    H = 0.5 * (H + H.transpose());

    VectorXd lower = VectorXd::Zero(variableDim);
    VectorXd upper = VectorXd::Zero(variableDim);
    const bool relaxLowerProgress =
        remaining <= _rmpcc.completionTolerance
        or remaining < static_cast<double>(N) * dt * progressRateMin
        or scheduleRemaining < static_cast<double>(N) * dt * progressRateMin;

    for(int stage = 0; stage < N; ++stage)
    {
        const int uOffset = stage * TWIST_DIM;
        lower.segment<3>(uOffset) = -_rmpcc.linearVelocityMax;
        upper.segment<3>(uOffset) = _rmpcc.linearVelocityMax;
        lower.segment<3>(uOffset + 3) = -_rmpcc.angularVelocityMax;
        upper.segment<3>(uOffset + 3) = _rmpcc.angularVelocityMax;
        lower(progressOffset + stage) = relaxLowerProgress ? 0.0 : progressRateMin;
        upper(progressOffset + stage) = progressRateMax;
    }

    MatrixXd Bineq = MatrixXd::Zero(2 * variableDim + 2, variableDim);
    VectorXd zineq = VectorXd::Zero(2 * variableDim + 2);
    Bineq.topRows(variableDim).setIdentity();
    zineq.head(variableDim) = upper;
    Bineq.block(variableDim, 0, variableDim, variableDim) = -MatrixXd::Identity(variableDim, variableDim);
    zineq.segment(variableDim, variableDim) = -lower;
    for(int stage = 0; stage < N; ++stage)
    {
        Bineq(2 * variableDim, progressOffset + stage) = dt;
        Bineq(2 * variableDim + 1, progressOffset + stage) = dt;
    }
    zineq(2 * variableDim) = remaining + _rmpcc.progressUpperSlack;
    zineq(2 * variableDim + 1) = scheduleRemaining;

    if(_warmStart.size() != variableDim)
    {
        _warmStart = VectorXd::Zero(variableDim);
        for(int stage = 0; stage < N; ++stage)
        {
            _warmStart.segment<TWIST_DIM>(stage * TWIST_DIM) =
                body_twist_reference_at_progress(predictedProgress[static_cast<size_t>(stage)]);
            _warmStart(progressOffset + stage) = _rmpcc.progressRateRef;
        }
    }

    VectorXd z0 = clipped_warm_start(_warmStart, lower, upper, dt, remaining, scheduleRemaining);
    VectorXd zOpt;
    bool fallbackUsed = false;
    double qpStatus = 1.0;
    try
    {
        zOpt = _qpSolver.solve(H, f, Bineq, zineq, z0);
        if(zOpt.size() != variableDim or not zOpt.allFinite())
        {
            throw std::runtime_error("QP returned an invalid solution.");
        }
    }
    catch(const std::exception &)
    {
        zOpt = z0;
        if(remaining <= _rmpcc.completionTolerance)
        {
            for(int stage = 0; stage < N; ++stage)
            {
                zOpt.segment<TWIST_DIM>(stage * TWIST_DIM).setZero();
            }
        }
        fallbackUsed = true;
        qpStatus = 0.0;
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
    _diagnostics.rotationError = rotationError;
    _diagnostics.qpStatus = qpStatus;
    _diagnostics.fallbackUsed = fallbackUsed;
    _diagnostics.referenceLinearSpeed = (currentTau * _rmpcc.progressRateRef).head<3>().norm();
    _diagnostics.referenceAngularSpeed = (currentTau * _rmpcc.progressRateRef).tail<3>().norm();
    const Eigen::Matrix<double, 1, TWIST_DIM> lagRow =
        (currentTau.transpose() * _rmpcc.metric) / std::sqrt(currentDenom);
    _diagnostics.lagError = (lagRow * e0)(0);
    _diagnostics.contourError = currentContourError;
}

Eigen::VectorXd
SerialLinkRMPCC::step(const double dt, const double estimatedProgress)
{
    if(not _trajectorySet)
    {
        throw std::runtime_error("[ERROR] [SERIAL LINK RMPCC] step(): No trajectory set. Call set_trajectory() first.");
    }
    if(not std::isfinite(dt) or dt <= 0.0)
    {
        throw std::invalid_argument("[ERROR] [SERIAL LINK RMPCC] step(): dt must be positive.");
    }
    if(not std::isfinite(estimatedProgress))
    {
        throw std::invalid_argument(
            "[ERROR] [SERIAL LINK RMPCC] step(): estimatedProgress must be finite.");
    }

    _pathProgress = clamp_value(estimatedProgress, 0.0, 1.0);
    update();

    const RobotLibrary::Model::Pose currentPose = endpoint_pose();
    const Eigen::Matrix4d currentTransform = pose_to_matrix(currentPose);

    solve_rmpcc(currentTransform, dt);

    // The model Jacobian maps qdot to endpoint point velocity [p_dot; omega],
    // not to the screw-theory spatial twist [v; omega]. Rotate both components
    // into the base frame without the SE(3) adjoint p x omega term.
    const Eigen::Matrix4d referenceTransform = reference_transform(_pathProgress);
    const Eigen::Matrix3d referenceRotation = referenceTransform.block<3,3>(0,0);
    Eigen::Vector<double, TWIST_DIM> baseTwist;
    baseTwist.head<3>() = referenceRotation * _diagnostics.bodyTwist.head<3>();
    baseTwist.tail<3>() = referenceRotation * _diagnostics.bodyTwist.tail<3>();

    if(_rmpcc.poseFeedbackEnable)
    {
        const RobotLibrary::Model::Pose referencePose = matrix_to_pose(referenceTransform);
        baseTwist += _rmpcc.poseFeedbackGain * pose_feedback_error(currentPose, referencePose);
    }

    _lastBodyTwist = _diagnostics.bodyTwist;
    _lastProgressRate = _diagnostics.progressRate;
    _diagnostics.pathProgress = _pathProgress;

    return resolve_endpoint_twist(baseTwist);
}

} } // namespace
