/**
 * @file    SerialLinkRMPCC.cpp
 * @brief   Riemannian MPCC controller implementation for serial link robot arms.
 */

#include <Control/TrajectoryTracking/SerialLinkRMPCC.h>
#include "detail/ProgressSchedule.h"
#include "detail/RmpccPrediction.h"
#include "detail/RmpccPhaseResidual.h"
#include "detail/RmpccAssociatedPhase.h"
#include "detail/RmpccProgressConstraints.h"
#include "detail/RmpccQpConstraints.h"
#include "detail/RmpccReferenceMotion.h"
#include "detail/RmpccResidualLinearization.h"
#include <Math/MathFunctions.h>

#include <Eigen/Geometry>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
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

template<typename Derived>
void hash_eigen(std::uint64_t &hash, const Eigen::MatrixBase<Derived> &value)
{
    constexpr std::uint64_t prime = 1099511628211ULL;
    const unsigned char *bytes = reinterpret_cast<const unsigned char *>(
        value.derived().data());
    const std::size_t byteCount = static_cast<std::size_t>(value.size())
        * sizeof(typename Derived::Scalar);
    for(std::size_t i = 0; i < byteCount; ++i)
    {
        hash ^= bytes[i];
        hash *= prime;
    }
}

RobotLibrary::Model::Pose matrix_to_pose(const Eigen::Matrix4d &T)
{
    Eigen::Quaterniond q(T.block<3,3>(0,0));
    q.normalize();
    return RobotLibrary::Model::Pose(T.block<3,1>(0,3), q);
}

template<typename Derived>
void write_csv(const std::filesystem::path &path,
               const Eigen::MatrixBase<Derived> &value)
{
    std::ofstream stream(path);
    if(not stream)
    {
        throw std::runtime_error("Unable to create N125-QP snapshot file: "
                                 + path.string());
    }
    stream << std::setprecision(17);
    for(Eigen::Index row = 0; row < value.rows(); ++row)
    {
        for(Eigen::Index column = 0; column < value.cols(); ++column)
        {
            if(column != 0)
            {
                stream << ',';
            }
            stream << value(row, column);
        }
        stream << '\n';
    }
}

std::string environment_value(const char *name)
{
    const char *value = std::getenv(name);
    return value == nullptr ? std::string() : std::string(value);
}

std::string json_escape(const std::string &value)
{
    std::ostringstream stream;
    for(const char character : value)
    {
        switch(character)
        {
            case '\\': stream << "\\\\"; break;
            case '"': stream << "\\\""; break;
            case '\n': stream << "\\n"; break;
            case '\r': stream << "\\r"; break;
            case '\t': stream << "\\t"; break;
            default: stream << character; break;
        }
    }
    return stream.str();
}

void write_n125_qp_snapshot(
    const Eigen::MatrixXd &H,
    const Eigen::VectorXd &f,
    const Eigen::MatrixXd &Aeq,
    const Eigen::VectorXd &yeq,
    const Eigen::MatrixXd &Bineq,
    const Eigen::VectorXd &zineq,
    const Eigen::VectorXd &zNominal,
    const Eigen::VectorXd &lower,
    const Eigen::VectorXd &upper,
    const Eigen::VectorXd &fixedProgressRates,
    const Eigen::Vector<double,6> &eta,
    const std::vector<RmpccStateVector> &nominalStates,
    const std::vector<RmpccInputVector> &nominalInputs,
    const Eigen::VectorXd &zReturned,
    const SolverResults<double> &solverResults,
    const std::uint64_t controlStepIndex,
    const double pathProgress,
    const double scheduleProgressLimit,
    const double remainingProgress,
    const double scheduleRemaining,
    const double controllerDt)
{
    const std::string requestedDirectory =
        environment_value("RMPCC_QP_SNAPSHOT_DIR");
    if(requestedDirectory.empty())
    {
        return;
    }

    const std::filesystem::path directory(requestedDirectory);
    std::filesystem::create_directories(directory);
    const std::filesystem::path metadataPath = directory / "metadata.json";
    if(std::filesystem::exists(metadataPath))
    {
        return;
    }

    Eigen::MatrixXd states(static_cast<Eigen::Index>(nominalStates.size()), 7);
    for(Eigen::Index row = 0; row < states.rows(); ++row)
    {
        states.row(row) = nominalStates[static_cast<std::size_t>(row)].transpose();
    }
    Eigen::MatrixXd inputs(static_cast<Eigen::Index>(nominalInputs.size()), 7);
    for(Eigen::Index row = 0; row < inputs.rows(); ++row)
    {
        inputs.row(row) = nominalInputs[static_cast<std::size_t>(row)].transpose();
    }

    write_csv(directory / "H.csv", H);
    write_csv(directory / "f.csv", f);
    write_csv(directory / "Aeq.csv", Aeq);
    write_csv(directory / "yeq.csv", yeq);
    write_csv(directory / "Bineq.csv", Bineq);
    write_csv(directory / "zineq.csv", zineq);
    write_csv(directory / "zNominal.csv", zNominal);
    write_csv(directory / "lower.csv", lower);
    write_csv(directory / "upper.csv", upper);
    write_csv(directory / "fixed_progress_rates.csv", fixedProgressRates);
    write_csv(directory / "eta.csv", eta);
    write_csv(directory / "nominal_states.csv", states);
    write_csv(directory / "nominal_inputs.csv", inputs);
    write_csv(directory / "zReturned.csv", zReturned);

    const double nominalViolation = (Bineq * zNominal - zineq).maxCoeff();
    const double returnedViolation = (Bineq * zReturned - zineq).maxCoeff();
    const auto timestamp = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::system_clock::now().time_since_epoch()).count();
    std::ofstream metadata(metadataPath);
    if(not metadata)
    {
        throw std::runtime_error("Unable to create N125-QP snapshot metadata: "
                                 + metadataPath.string());
    }
    metadata << std::setprecision(17)
             << "{\n"
             << "  \"schema\": \"n125_qp_snapshot_v2\",\n"
             << "  \"timestamp_unix_ns\": " << timestamp << ",\n"
             << "  \"control_step_index_zero_based\": " << controlStepIndex << ",\n"
             << "  \"profile\": \""
             << json_escape(environment_value("RMPCC_QP_PROFILE")) << "\",\n"
             << "  \"controlled_frame\": \""
             << json_escape(environment_value("RMPCC_QP_CONTROLLED_FRAME")) << "\",\n"
             << "  \"trajectory\": \""
             << json_escape(environment_value("RMPCC_QP_TRAJECTORY")) << "\",\n"
             << "  \"controller_dt_s\": " << controllerDt << ",\n"
             << "  \"horizon_steps\": " << fixedProgressRates.size() << ",\n"
             << "  \"variable_dim\": " << zNominal.size() << ",\n"
             << "  \"equality_rows\": " << Aeq.rows() << ",\n"
             << "  \"inequality_rows\": " << Bineq.rows() << ",\n"
             << "  \"current_path_progress\": " << pathProgress << ",\n"
             << "  \"schedule_progress_limit\": " << scheduleProgressLimit << ",\n"
             << "  \"remaining_progress\": " << remainingProgress << ",\n"
             << "  \"schedule_remaining\": " << scheduleRemaining << ",\n"
             << "  \"nominal_max_violation\": " << nominalViolation << ",\n"
             << "  \"returned_max_violation\": " << returnedViolation << ",\n"
             << "  \"solver_number_of_steps\": " << solverResults.numberOfSteps << ",\n"
             << "  \"solver_final_step_size\": " << solverResults.finalStepSize << ",\n"
             << "  \"solver_objective\": " << solverResults.objectiveFunction << "\n"
             << "}\n";
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
  _qpSolver(parameters.qpsolver),
  _poseArcTable(std::make_shared<RmpccPoseArcTable>())
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
    if(_rmpcc.referenceMotion == RmpccReferenceMotion::StageConsistent
       && (not std::isfinite(_rmpcc.geometricTangentStep)
           or _rmpcc.geometricTangentStep <= 0.0))
    {
        throw std::invalid_argument(
            "[ERROR] [SERIAL LINK RMPCC] StageConsistent reference motion requires a positive geometricTangentStep.");
    }
    if(_rmpcc.referenceMotion == RmpccReferenceMotion::StageConsistent
       && _rmpcc.predictorGeometry != RmpccPredictorGeometry::ExactSE3)
    {
        throw std::invalid_argument(
            "[ERROR] [SERIAL LINK RMPCC] StageConsistent reference motion requires the ExactSE3 predictor.");
    }
    if(not std::isfinite(_rmpcc.phaseDenominatorTolerance)
       or _rmpcc.phaseDenominatorTolerance < 0.0)
    {
        throw std::invalid_argument(
            "[ERROR] [SERIAL LINK RMPCC] phaseDenominatorTolerance must be finite and non-negative.");
    }
    const bool associatedPhase =
        _rmpcc.contourResidualGeometry
            == RmpccContourResidualGeometry::ExactAssociatedUnifiedSE3
        or _rmpcc.contourResidualGeometry
            == RmpccContourResidualGeometry::AssociatedDecoupledCartesianSO3;
    if((_rmpcc.phaseAssociation == RmpccPhaseAssociation::TaskPointXYZ
        or _rmpcc.phaseAssociation == RmpccPhaseAssociation::TaskPoseFeature
        or _rmpcc.lagPenalty == RmpccLagPenalty::ScalarTaskDistance
        or _rmpcc.lagPenalty == RmpccLagPenalty::ScalarPosePathArc
        or associatedPhase)
       and _rmpcc.residualLinearization
               != RmpccResidualLinearization::FullResidualJacobian)
    {
        throw std::invalid_argument(
            "[ERROR] [SERIAL LINK RMPCC] Task-associated/scalar lag modes require FullResidualJacobian.");
    }
    if(_rmpcc.phaseAssociation == RmpccPhaseAssociation::TaskPoseFeature
       && (not std::isfinite(_rmpcc.rotationCharacteristicLength)
           or _rmpcc.rotationCharacteristicLength <= 0.0))
    {
        throw std::invalid_argument(
            "[ERROR] [SERIAL LINK RMPCC] TaskPoseFeature requires a positive rotationCharacteristicLength.");
    }
    if(_rmpcc.lagPenalty == RmpccLagPenalty::ScalarTaskDistance)
    {
        const Eigen::Vector3d translationWeights =
            _rmpcc.lagWeightMatrix.diagonal().head<3>();
        if((translationWeights.array() - translationWeights(0)).abs().maxCoeff()
           > 1e-12 * std::max(1.0, std::abs(translationWeights(0))))
        {
            throw std::invalid_argument(
                "[ERROR] [SERIAL LINK RMPCC] ScalarTaskDistance requires equal translational lag weights.");
        }
    }
    if(associatedPhase
       && _rmpcc.phaseAssociation != RmpccPhaseAssociation::TaskPointXYZ
       && _rmpcc.phaseAssociation != RmpccPhaseAssociation::TaskPoseFeature)
    {
        throw std::invalid_argument(
            "[ERROR] [SERIAL LINK RMPCC] Associated contour residuals require a task phase association.");
    }
    // TaskPoseFeature only reaches the cost through the associated family or
    // through the full-screw phase family. Anywhere else it is silently inert,
    // which would make an experiment look configured when it is not.
    if(_rmpcc.phaseAssociation == RmpccPhaseAssociation::TaskPoseFeature
       && not associatedPhase
       && (_rmpcc.objectiveGeometry != RmpccObjectiveGeometry::FullScrewSE3
           or _rmpcc.runningLagGeometry != RmpccLagGeometry::FullScrew))
    {
        throw std::invalid_argument(
            "[ERROR] [SERIAL LINK RMPCC] TaskPoseFeature has no cost effect for this objective "
            "(it requires FullScrewSE3 with FullScrew running lag, or an associated contour residual).");
    }
    if(_rmpcc.contourResidualGeometry
           == RmpccContourResidualGeometry::AssociatedDecoupledCartesianSO3
       && (_rmpcc.objectiveGeometry
               != RmpccObjectiveGeometry::DecoupledCartesianSO3
           or _rmpcc.lagPenalty != RmpccLagPenalty::ScalarPosePathArc))
    {
        throw std::invalid_argument(
            "[ERROR] [SERIAL LINK RMPCC] Associated decoupled mode requires DecoupledCartesianSO3 and ScalarPosePathArc.");
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
    validateCostScale(_rmpcc.trackingCostScale, "trackingCostScale");
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
        stream << "phase=";
        if(_rmpcc.phaseAssociation == RmpccPhaseAssociation::MetricScrew)
            stream << "metric_screw";
        else if(_rmpcc.phaseAssociation == RmpccPhaseAssociation::TaskPointXYZ)
            stream << "task_point_xyz";
        else
            stream << "task_pose_feature";
        stream
               << "; lag_penalty="
               << (_rmpcc.lagPenalty == RmpccLagPenalty::PhaseInducedPoseVector
                   ? "phase_induced_pose_vector"
                   : (_rmpcc.lagPenalty == RmpccLagPenalty::ScalarTaskDistance
                      ? "scalar_task_distance" : "scalar_pose_path_arc"))
               << "; contour="
               << (_rmpcc.contourResidualGeometry
                       == RmpccContourResidualGeometry::ExactAssociatedUnifiedSE3
                   ? "Log(Tref(s_hat)^-1*Tactual)" : "e-g*delta_s");
    }
    else
    {
        stream << "sum ||(I-tp*tp^T)e_p||_Qcp^2+||tp*tp^T*e_p||_Qlp^2"
               << "+||Log(Rref^T*R)||_QR^2; orientation is not projected onto progress";
    }
    if(_rmpcc.contourResidualGeometry
           == RmpccContourResidualGeometry::AssociatedDecoupledCartesianSO3)
    {
        stream << "; associated_phase="
               << (_rmpcc.phaseAssociation == RmpccPhaseAssociation::TaskPoseFeature
                   ? "task_pose_feature" : "task_point_xyz")
               << "; contour_reference=s_hat"
               << "; lag_penalty=scalar_pose_path_arc";
    }
    stream << "; reference_motion="
           << (_rmpcc.referenceMotion == RmpccReferenceMotion::StageConsistent
               ? "stage_consistent" : "legacy_tangent_product")
           << "; path_velocity="
           << (_rmpcc.referenceMotion == RmpccReferenceMotion::StageConsistent
               ? "sum ||u-Ad(E^-1)*Log(T(s)^-1*T(s+dt*sdot))/dt||_Rv^2"
               : "sum ||u-Ad(E^-1)*tau*sdot||_Rv^2")
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

    const bool associatedPhase =
        _rmpcc.contourResidualGeometry
            == RmpccContourResidualGeometry::ExactAssociatedUnifiedSE3
        or _rmpcc.contourResidualGeometry
            == RmpccContourResidualGeometry::AssociatedDecoupledCartesianSO3;
    if(_rmpcc.lagPenalty == RmpccLagPenalty::ScalarPosePathArc
       or associatedPhase)
    {
        _poseArcTable->build(
            _rmpcc.lagWeightMatrix,
            [this](const double progress)
            {
                if(_rmpcc.referenceMotion
                   == RmpccReferenceMotion::LegacyTangentProduct)
                {
                    return _trajectory.tangent_at_progress(
                        progress, _rmpcc.tangentStep);
                }
                const auto reference = [this](const double s)
                {
                    return _trajectory.pose_at_progress(s).as_matrix();
                };
                return rmpcc_centred_geometric_tangent(
                    progress, _rmpcc.geometricTangentStep, reference);
            },
            4097);
    }

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
    _parentFrameMotion.set_static_pose(active_frame_transform_in_base());
}

void
SerialLinkRMPCC::set_trajectory_frame(
    const RobotLibrary::Trajectory::CartesianTrajectoryFrameState &frame,
    const double timestampSeconds)
{
    RobotLibrary::Trajectory::validate_trajectory_frame(frame);
    _trajectoryFrame = frame;
    _parentFrameMotion.update(active_frame_transform_in_base(), timestampSeconds);
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
    _controlStepIndex = 0;
    _scheduleProgressLimit = 1.0;
    _parentFrameMotion.reset();
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

    const auto totalSolveStart = std::chrono::steady_clock::now();
    double referencePoseQueryTime = 0.0;
    double referenceTangentQueryTime = 0.0;
    std::uint64_t referencePoseQueryCount = 0;
    std::uint64_t referenceTangentQueryCount = 0;
    std::uint64_t referencePoseRequestCount = 0;
    std::uint64_t referenceTangentRequestCount = 0;
    double rolloutLinearizationTime = 0.0;
    double residualHessianTime = 0.0;
    const bool numericalAudit = std::getenv("ROBOT_LIBRARY_RMPCC_NUMERICAL_AUDIT") != nullptr;
    std::uint64_t stateLinearizationHash = 1469598103934665603ULL;
    std::unordered_map<std::uint64_t, Eigen::Matrix4d> referenceTransformCache;
    std::unordered_map<std::uint64_t, Eigen::Vector<double,TWIST_DIM>> referenceTangentCache;
    const auto progressKey = [](const double progress)
    {
        std::uint64_t key = 0;
        static_assert(sizeof(key) == sizeof(progress));
        std::memcpy(&key, &progress, sizeof(key));
        return key;
    };
    const auto referenceTransform = [&](const double progress)
    {
        const auto start = std::chrono::steady_clock::now();
        ++referencePoseRequestCount;
        const std::uint64_t key = progressKey(progress);
        const auto cached = referenceTransformCache.find(key);
        if(cached != referenceTransformCache.end())
        {
            referencePoseQueryTime += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start).count();
            return cached->second;
        }
        const Eigen::Matrix4d value =
            reference_transform_in_trajectory_frame(progress);
        referenceTransformCache.emplace(key, value);
        referencePoseQueryTime += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        ++referencePoseQueryCount;
        return value;
    };
    const auto referenceTangent = [&](const double progress)
    {
        const auto start = std::chrono::steady_clock::now();
        ++referenceTangentRequestCount;
        const std::uint64_t key = progressKey(progress);
        const auto cached = referenceTangentCache.find(key);
        if(cached != referenceTangentCache.end())
        {
            referenceTangentQueryTime += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start).count();
            return cached->second;
        }
        const Eigen::Vector<double,TWIST_DIM> value =
            _rmpcc.referenceMotion == RmpccReferenceMotion::LegacyTangentProduct
            ? _trajectory.tangent_at_progress(progress, _rmpcc.tangentStep)
            : rmpcc_centred_geometric_tangent(
                  progress, _rmpcc.geometricTangentStep, referenceTransform);
        referenceTangentCache.emplace(key, value);
        referenceTangentQueryTime += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        ++referenceTangentQueryCount;
        return value;
    };
    const bool parentMotionActive = _parentFrameMotion.has_velocity()
        && _parentFrameMotion.body_twist().squaredNorm() > 0.0;
    const bool parentRepairActive = parentMotionActive
        && _rmpcc.referenceMotion == RmpccReferenceMotion::LegacyTangentProduct;
    const auto parentTransform = [&](const int stage)
    {
        return _parentFrameMotion.predicted_pose(stage, dt);
    };

    _diagnostics.qpStatus = 0.0;
    _diagnostics.fallbackUsed = false;
    const std::uint64_t controlStepIndex = _controlStepIndex++;

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
        referenceTransform(_pathProgress);
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
        referenceTangent(_pathProgress);
    const Eigen::Vector<double, TWIST_DIM> currentErrorTangent =
        rmpcc_error_coordinate_path_tangent(e0, currentTau);
    RmpccStateVector currentState = RmpccStateVector::Zero();
    currentState.head<TWIST_DIM>() = e0;
    currentState(6) = _pathProgress;
    const auto currentReferenceTransformFunction = [&](const double progress)
    {
        return referenceTransform(progress);
    };
    const auto currentReferenceTangentFunction = [&](const double progress)
    {
        return referenceTangent(progress);
    };
    const RmpccPhaseResiduals currentPhaseResidual = rmpcc_phase_residuals(
        currentState, _rmpcc.metric, _rmpcc.phaseAssociation,
        _rmpcc.hessianRegularization, _rmpcc.phaseDenominatorTolerance,
        _rmpcc.rotationCharacteristicLength,
        currentReferenceTransformFunction, currentReferenceTangentFunction);
    const bool associatedPhaseFamily =
        _rmpcc.contourResidualGeometry
            == RmpccContourResidualGeometry::ExactAssociatedUnifiedSE3
        or _rmpcc.contourResidualGeometry
            == RmpccContourResidualGeometry::AssociatedDecoupledCartesianSO3;
    RmpccAssociatedResiduals currentAssociatedResidual;
    if(associatedPhaseFamily)
    {
        currentAssociatedResidual = rmpcc_associated_residuals(
            currentState, _rmpcc.contourResidualGeometry, _rmpcc.metric,
            _rmpcc.phaseAssociation, _rmpcc.hessianRegularization,
            _rmpcc.phaseDenominatorTolerance,
            _rmpcc.rotationCharacteristicLength, *_poseArcTable,
            currentReferenceTransformFunction, currentReferenceTangentFunction);
    }
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
    const bool phaseFactorialFamily =
        _rmpcc.objectiveGeometry == RmpccObjectiveGeometry::FullScrewSE3
        && _rmpcc.runningLagGeometry == RmpccLagGeometry::FullScrew
        && not associatedPhaseFamily;
    const double currentContourError = associatedPhaseFamily
        ? currentAssociatedResidual.contour.norm()
        : (phaseFactorialFamily ? currentPhaseResidual.contour.norm()
                                : (currentContourProjection * currentCostError).norm());
    const double progressRateMax = _rmpcc.progressRateMax;
    const double progressRateMin = std::min(_rmpcc.progressRateMin, progressRateMax);

    const auto constraintStart = std::chrono::steady_clock::now();
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

    const RmpccQpConstraints qpConstraints = rmpcc_build_qp_constraints(
        lower, upper, N, dt, remaining + _rmpcc.progressUpperSlack,
        scheduleRemaining, _fixedProgressSchedule, fixedProgressRates);
    const MatrixXd &Aeq = qpConstraints.Aeq;
    const VectorXd &yeq = qpConstraints.yeq;
    const MatrixXd &Bineq = qpConstraints.Bineq;
    const VectorXd &zineq = qpConstraints.zineq;
    _diagnostics.constraintConstructionTimeSeconds =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - constraintStart).count();
    const auto propagate = [&](const RmpccStateVector &state,
                               const RmpccInputVector &input,
                               const int stage)
    {
        if(parentRepairActive)
        {
            return rmpcc_parent_repaired_legacy_state_step(
                state, input, dt, stage,
                referenceTransform, referenceTangent, parentTransform);
        }
        return _rmpcc.predictorGeometry == RmpccPredictorGeometry::ExactSE3
            ? rmpcc_exact_state_step(state, input, dt, referenceTransform)
            : rmpcc_additive_state_step(state, input, dt, referenceTangent);
    };
    const auto linearize = [&](const RmpccStateVector &state,
                               const RmpccInputVector &input,
                               const int stage)
    {
        if(parentRepairActive)
        {
            return rmpcc_linearize_parent_repaired_legacy_state_step(
                state, input, dt, stage, _rmpcc.rtiFiniteDifferenceStep,
                referenceTransform, referenceTangent, parentTransform);
        }
        return _rmpcc.predictorGeometry == RmpccPredictorGeometry::ExactSE3
            ? rmpcc_linearize_exact_state_step(
                  state, input, dt, _rmpcc.rtiFiniteDifferenceStep,
                  referenceTransform)
            : rmpcc_linearize_additive_state_step(
                  state, input, dt, _rmpcc.rtiFiniteDifferenceStep,
                  referenceTangent);
    };

    _diagnostics.warmStartInvariantError = 0.0;
    _diagnostics.warmStartInputBoundActive = 0.0;
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
            Eigen::Vector<double,TWIST_DIM> bodyTwist;
            if(parentRepairActive)
            {
                bodyTwist = rmpcc_parent_repaired_legacy_feedforward(
                    guessState.head<TWIST_DIM>(), guessState(6), rate, dt, stage,
                    referenceTransform, referenceTangent, parentTransform);
            }
            else if(_rmpcc.referenceMotion
               == RmpccReferenceMotion::StageConsistent)
            {
                bodyTwist = rmpcc_stage_consistent_feedforward(
                    guessState.head<TWIST_DIM>(), guessState(6), rate, dt,
                    referenceTransform);
            }
            else
            {
                const Eigen::Vector<double,TWIST_DIM> tangent =
                    referenceTangent(guessState(6));
                bodyTwist = rmpcc_transport_reference_tangent(
                    guessState.head<TWIST_DIM>(), tangent) * rate;
            }
            for(int component = 0; component < TWIST_DIM; ++component)
            {
                if(bodyTwist(component) < lower(uOffset + component)
                   or bodyTwist(component) > upper(uOffset + component))
                {
                    _diagnostics.warmStartInputBoundActive = 1.0;
                }
                bodyTwist(component) = clamp_value(bodyTwist(component),
                                                   lower(uOffset + component),
                                                   upper(uOffset + component));
            }
            _warmStart.segment<TWIST_DIM>(uOffset) = bodyTwist;
            _warmStart(sOffset) = rate;

            RmpccInputVector guessInput = RmpccInputVector::Zero();
            guessInput.head<TWIST_DIM>() = bodyTwist;
            guessInput(6) = rate;
            const Eigen::Matrix4d previousRelative =
                RobotLibrary::Math::se3_exponential(guessState.head<TWIST_DIM>());
            const RmpccStateVector nextGuessState =
                propagate(guessState, guessInput, stage);
            const Eigen::Matrix4d nextRelative =
                RobotLibrary::Math::se3_exponential(nextGuessState.head<TWIST_DIM>());
            _diagnostics.warmStartInvariantError = std::max(
                _diagnostics.warmStartInvariantError,
                RobotLibrary::Math::se3_logarithm(
                    RobotLibrary::Math::se3_inverse(previousRelative)
                    * nextRelative).norm());
            guessState = nextGuessState;
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
    MatrixXd scalarLagResidualSensitivity = MatrixXd::Zero(N, variableDim);
    VectorXd scalarLagResidualOffset = VectorXd::Zero(N);
    MatrixXd pathVelocityResidualSensitivity = MatrixXd::Zero(errorDim, variableDim);
    VectorXd pathVelocityResidualOffset = VectorXd::Zero(errorDim);
    Eigen::Matrix<double,7,Eigen::Dynamic> stateSensitivity(7, variableDim);
    stateSensitivity.setZero();

    std::vector<RmpccStateVector> nominalStates(static_cast<size_t>(N + 1));
    std::vector<RmpccInputVector> nominalInputs(static_cast<size_t>(N));
    std::vector<MatrixXd> stageStartStateSensitivities(
        static_cast<size_t>(N), MatrixXd::Zero(7, variableDim));
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
    std::vector<double> scalarLagCostWeights(static_cast<size_t>(N), 0.0);
    std::vector<bool> usesScalarLag(static_cast<size_t>(N), false);
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
        stageStartStateSensitivities[static_cast<size_t>(stage)] = stateSensitivity;

        const auto rolloutStart = std::chrono::steady_clock::now();
        const RmpccStageLinearization linearization = linearize(
            nominalStates[static_cast<size_t>(stage)], nominalInput, stage);
        if(numericalAudit)
        {
            hash_eigen(stateLinearizationHash, linearization.nominalNext);
            hash_eigen(stateLinearizationHash, linearization.stateJacobian);
            hash_eigen(stateLinearizationHash, linearization.inputJacobian);
        }
        nominalStates[static_cast<size_t>(stage + 1)] = linearization.nominalNext;

        stateSensitivity = linearization.stateJacobian * stateSensitivity;
        for(int component = 0; component < TWIST_DIM; ++component)
        {
            stateSensitivity.col(uOffset + component) +=
                linearization.inputJacobian.col(component);
        }
        stateSensitivity.col(sOffset) += linearization.inputJacobian.col(6);
        rolloutLinearizationTime += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - rolloutStart).count();
        const auto residualStart = std::chrono::steady_clock::now();

        const int row = stage * TWIST_DIM;
        const MatrixXd stageSensitivity = stateSensitivity.topRows(TWIST_DIM);

        const double predictedProgress = linearization.nominalNext(6);
        const Eigen::Vector<double,TWIST_DIM> stageReferenceTangent =
            referenceTangent(predictedProgress);
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
        stageContourWeight *= _rmpcc.trackingCostScale;
        stageLagWeight *= _rmpcc.trackingCostScale;
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
        scalarLagCostWeights[static_cast<size_t>(stage)] = stageLagWeight(0,0);
        const bool useFullResidualJacobian =
            _rmpcc.residualLinearization
                == RmpccResidualLinearization::FullResidualJacobian;
        if(useFullResidualJacobian)
        {
            MatrixXd contourJacobian;
            Eigen::Vector<double,TWIST_DIM> contourOffset;
            if(associatedPhaseFamily)
            {
                const RmpccAssociatedResidualLinearization residualLinearization =
                    rmpcc_linearize_associated_residuals(
                        linearization.nominalNext,
                        _rmpcc.contourResidualGeometry, _rmpcc.metric,
                        _rmpcc.phaseAssociation, _rmpcc.hessianRegularization,
                        _rmpcc.phaseDenominatorTolerance,
                        _rmpcc.rotationCharacteristicLength,
                        _rmpcc.rtiFiniteDifferenceStep, *_poseArcTable,
                        referenceTransform, referenceTangent);
                contourJacobian =
                    residualLinearization.contourJacobian * stateSensitivity;
                contourOffset = residualLinearization.residual.contour
                    - contourJacobian * zNominal;
                fullResidualH += 2.0 * contourJacobian.transpose()
                    * stageContourWeight * contourJacobian;
                fullResidualF += 2.0 * contourJacobian.transpose()
                    * stageContourWeight * contourOffset;
                if(_rmpcc.lagPenalty
                   == RmpccLagPenalty::PhaseInducedPoseVector)
                {
                    const MatrixXd lagJacobian =
                        residualLinearization.vectorLagJacobian * stateSensitivity;
                    const Eigen::Vector<double,TWIST_DIM> lagOffset =
                        residualLinearization.residual.vectorLag
                        - lagJacobian * zNominal;
                    fullResidualH += 2.0 * lagJacobian.transpose()
                        * stageLagWeight * lagJacobian;
                    fullResidualF += 2.0 * lagJacobian.transpose()
                        * stageLagWeight * lagOffset;
                    lagResidualSensitivity.block(
                        row, 0, TWIST_DIM, variableDim) = lagJacobian;
                    lagResidualOffset.segment<TWIST_DIM>(row) = lagOffset;
                }
                else
                {
                    const Eigen::Matrix<double,1,Eigen::Dynamic> scalarJacobian =
                        residualLinearization.scalarPoseArcLagJacobian
                        * stateSensitivity;
                    const double scalarOffset =
                        residualLinearization.residual.scalarPoseArcLag
                        - (scalarJacobian * zNominal)(0);
                    const double qLag =
                        scalarLagCostWeights[static_cast<size_t>(stage)];
                    fullResidualH +=
                        2.0 * qLag * scalarJacobian.transpose() * scalarJacobian;
                    fullResidualF +=
                        2.0 * qLag * scalarJacobian.transpose() * scalarOffset;
                    scalarLagResidualSensitivity.row(stage) = scalarJacobian;
                    scalarLagResidualOffset(stage) = scalarOffset;
                    usesScalarLag[static_cast<size_t>(stage)] = true;
                }
            }
            else if(phaseFactorialFamily)
            {
                const RmpccPhaseResidualLinearization residualLinearization =
                    rmpcc_linearize_phase_residuals(
                        linearization.nominalNext, _rmpcc.metric,
                        _rmpcc.phaseAssociation, _rmpcc.hessianRegularization,
                        _rmpcc.phaseDenominatorTolerance,
                        _rmpcc.rotationCharacteristicLength,
                        _rmpcc.rtiFiniteDifferenceStep,
                        referenceTransform, referenceTangent);
                contourJacobian =
                    residualLinearization.contourJacobian * stateSensitivity;
                contourOffset = residualLinearization.residual.contour
                    - contourJacobian * zNominal;
                fullResidualH += 2.0 * contourJacobian.transpose()
                    * stageContourWeight * contourJacobian;
                fullResidualF += 2.0 * contourJacobian.transpose()
                    * stageContourWeight * contourOffset;
                if(_rmpcc.lagPenalty
                   == RmpccLagPenalty::PhaseInducedPoseVector)
                {
                    const MatrixXd lagJacobian =
                        residualLinearization.vectorLagJacobian * stateSensitivity;
                    const Eigen::Vector<double,TWIST_DIM> lagOffset =
                        residualLinearization.residual.vectorLag
                        - lagJacobian * zNominal;
                    fullResidualH += 2.0 * lagJacobian.transpose()
                        * stageLagWeight * lagJacobian;
                    fullResidualF += 2.0 * lagJacobian.transpose()
                        * stageLagWeight * lagOffset;
                    lagResidualSensitivity.block(
                        row, 0, TWIST_DIM, variableDim) = lagJacobian;
                    lagResidualOffset.segment<TWIST_DIM>(row) = lagOffset;
                }
                else
                {
                    const Eigen::Matrix<double,1,Eigen::Dynamic> scalarJacobian =
                        residualLinearization.scalarLagJacobian * stateSensitivity;
                    const double scalarOffset =
                        residualLinearization.residual.scalarLag
                        - (scalarJacobian * zNominal)(0);
                    const double qLag = scalarLagCostWeights[static_cast<size_t>(stage)];
                    fullResidualH +=
                        2.0 * qLag * scalarJacobian.transpose() * scalarJacobian;
                    fullResidualF +=
                        2.0 * qLag * scalarJacobian.transpose() * scalarOffset;
                    scalarLagResidualSensitivity.row(stage) = scalarJacobian;
                    scalarLagResidualOffset(stage) = scalarOffset;
                    usesScalarLag[static_cast<size_t>(stage)] = true;
                }
            }
            else
            {
                const RmpccFullScrewResidualLinearization residualLinearization =
                    _rmpcc.objectiveGeometry == RmpccObjectiveGeometry::FullScrewSE3
                    ? rmpcc_linearize_projected_residuals(
                          linearization.nominalNext, _rmpcc.metric, lagGeometry,
                          _rmpcc.hessianRegularization,
                          _rmpcc.rtiFiniteDifferenceStep, referenceTangent)
                    : rmpcc_linearize_decoupled_residuals(
                          linearization.nominalNext, _rmpcc.hessianRegularization,
                          _rmpcc.rtiFiniteDifferenceStep, referenceTangent);
                contourJacobian =
                    residualLinearization.contourJacobian * stateSensitivity;
                const MatrixXd lagJacobian =
                    residualLinearization.lagJacobian * stateSensitivity;
                contourOffset = residualLinearization.residual.contour
                    - contourJacobian * zNominal;
                const Eigen::Vector<double,TWIST_DIM> lagOffset =
                    residualLinearization.residual.lag - lagJacobian * zNominal;
                fullResidualH += 2.0 * contourJacobian.transpose()
                    * stageContourWeight * contourJacobian
                    + 2.0 * lagJacobian.transpose()
                    * stageLagWeight * lagJacobian;
                fullResidualF += 2.0 * contourJacobian.transpose()
                    * stageContourWeight * contourOffset
                    + 2.0 * lagJacobian.transpose()
                    * stageLagWeight * lagOffset;
                lagResidualSensitivity.block(row, 0, TWIST_DIM, variableDim) =
                    lagJacobian;
                lagResidualOffset.segment<TWIST_DIM>(row) = lagOffset;
            }
            contourResidualSensitivity.block(row, 0, TWIST_DIM, variableDim) =
                contourJacobian;
            contourResidualOffset.segment<TWIST_DIM>(row) = contourOffset;
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
        residualHessianTime += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - residualStart).count();
    }

    const auto finalAssemblyStart = std::chrono::steady_clock::now();
    MatrixXd H = 2.0 * errorSensitivity.transpose() * Q * errorSensitivity
                 + fullResidualH;
    VectorXd f = 2.0 * errorSensitivity.transpose() * Q * errorOffset
                 + fullResidualF;

    for(int stage = 0; stage < N; ++stage)
    {
        const int uOffset = stage * TWIST_DIM;
        const int sOffset = progressOffset + stage;
        const RmpccStateVector &stageState = nominalStates[static_cast<size_t>(stage)];
        const Eigen::Matrix<double,TWIST_DIM,TWIST_DIM> pathVelocityWeight =
            _rmpcc.pathVelocityScale * _rmpcc.pathVelocityWeight;
        if(_rmpcc.referenceMotion == RmpccReferenceMotion::StageConsistent
           || parentRepairActive)
        {
            H.block<TWIST_DIM,TWIST_DIM>(uOffset, uOffset) +=
                2.0 * _rmpcc.controlWeight;
            const RmpccPathVelocityResidualLinearization residualLinearization =
                parentRepairActive
                ? rmpcc_linearize_parent_repaired_legacy_path_velocity_residual(
                      stageState, nominalInputs[static_cast<size_t>(stage)], dt,
                      stage, _rmpcc.rtiFiniteDifferenceStep,
                      referenceTransform, referenceTangent, parentTransform)
                : rmpcc_linearize_stage_consistent_path_velocity_residual(
                      stageState, nominalInputs[static_cast<size_t>(stage)], dt,
                      _rmpcc.rtiFiniteDifferenceStep, referenceTransform);
            MatrixXd residualJacobian = residualLinearization.stateJacobian
                * stageStartStateSensitivities[static_cast<size_t>(stage)];
            residualJacobian.block<TWIST_DIM,TWIST_DIM>(0, uOffset) +=
                residualLinearization.inputJacobian.leftCols<TWIST_DIM>();
            residualJacobian.col(sOffset) +=
                residualLinearization.inputJacobian.col(6);
            const Eigen::Vector<double,TWIST_DIM> residualOffset =
                residualLinearization.residual - residualJacobian * zNominal;
            H += 2.0 * residualJacobian.transpose()
                * pathVelocityWeight * residualJacobian;
            f += 2.0 * residualJacobian.transpose()
                * pathVelocityWeight * residualOffset;
            pathVelocityResidualSensitivity.block(
                stage * TWIST_DIM, 0, TWIST_DIM, variableDim) = residualJacobian;
            pathVelocityResidualOffset.segment<TWIST_DIM>(stage * TWIST_DIM) =
                residualOffset;
        }
        else
        {
            const Eigen::Vector<double,TWIST_DIM> tangent =
                referenceTangent(stageState(6));
            const Eigen::Vector<double,TWIST_DIM> transportedTangent =
                rmpcc_transport_reference_tangent(
                    stageState.head<TWIST_DIM>(), tangent);
            transportedTangents[static_cast<size_t>(stage)] = transportedTangent;
            H.block<TWIST_DIM,TWIST_DIM>(uOffset, uOffset) +=
                2.0 * (_rmpcc.controlWeight + pathVelocityWeight);
            H.block<TWIST_DIM,1>(uOffset, sOffset) +=
                -2.0 * pathVelocityWeight * transportedTangent;
            H.block<1,TWIST_DIM>(sOffset, uOffset) +=
                -2.0 * transportedTangent.transpose() * pathVelocityWeight;
            H(sOffset, sOffset) +=
                2.0 * (transportedTangent.transpose()
                       * pathVelocityWeight * transportedTangent)(0);
        }

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
    if(numericalAudit)
    {
        std::uint64_t residualHash = 1469598103934665603ULL;
        hash_eigen(residualHash, errorSensitivity);
        hash_eigen(residualHash, errorOffset);
        hash_eigen(residualHash, contourResidualSensitivity);
        hash_eigen(residualHash, contourResidualOffset);
        hash_eigen(residualHash, lagResidualSensitivity);
        hash_eigen(residualHash, lagResidualOffset);
        hash_eigen(residualHash, scalarLagResidualSensitivity);
        hash_eigen(residualHash, scalarLagResidualOffset);
        if(_rmpcc.referenceMotion == RmpccReferenceMotion::StageConsistent)
        {
            hash_eigen(residualHash, pathVelocityResidualSensitivity);
            hash_eigen(residualHash, pathVelocityResidualOffset);
        }
        _diagnostics.stateLinearizationHash = stateLinearizationHash;
        _diagnostics.residualLinearizationHash = residualHash;
        std::uint64_t hessianHash = 1469598103934665603ULL;
        hash_eigen(hessianHash, H);
        hash_eigen(hessianHash, f);
        _diagnostics.hessianHash = hessianHash;
    }
    residualHessianTime += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - finalAssemblyStart).count();

    const auto solveStart = std::chrono::steady_clock::now();
    VectorXd zOpt = _fixedProgressSchedule
        ? _qpSolver.solve(H, f, Aeq, yeq, Bineq, zineq, zNominal)
        : _qpSolver.solve(H, f, Bineq, zineq, zNominal);
    _diagnostics.qpSolveTimeSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - solveStart).count();
    if(zOpt.size() != variableDim or not zOpt.allFinite())
    {
        throw std::runtime_error("[ERROR] [SERIAL LINK RMPCC] solve_rmpcc(): QP returned an invalid solution.");
    }
    const double inequalityViolation = (Bineq * zOpt - zineq).maxCoeff();
    const double equalityViolation = Aeq.rows() == 0
        ? 0.0 : (Aeq * zOpt - yeq).cwiseAbs().maxCoeff();
    const double constraintViolation =
        std::max(inequalityViolation, equalityViolation);
    _diagnostics.qpPrimalViolation = constraintViolation;
    _diagnostics.qpEqualityRows = static_cast<double>(Aeq.rows());
    _diagnostics.qpInequalityRows = static_cast<double>(Bineq.rows());
    _diagnostics.qpEqualityResidual = equalityViolation;
    _diagnostics.qpInequalityViolation = std::max(0.0, inequalityViolation);
    if(not std::isfinite(constraintViolation) or constraintViolation > 1e-6)
    {
        write_n125_qp_snapshot(
            H, f, Aeq, yeq, Bineq, zineq, zNominal, lower, upper, fixedProgressRates,
            e0, nominalStates, nominalInputs, zOpt, _qpSolver.results(),
            controlStepIndex, _pathProgress, _scheduleProgressLimit, remaining,
            scheduleRemaining, dt);
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
        parentRepairActive
        ? rmpcc_parent_repaired_legacy_feedforward(
              e0, _pathProgress, _diagnostics.progressRate, dt, 0,
              referenceTransform, referenceTangent, parentTransform)
        : _rmpcc.referenceMotion == RmpccReferenceMotion::StageConsistent
        ? rmpcc_stage_consistent_feedforward(
              e0, _pathProgress, _diagnostics.progressRate, dt,
              referenceTransform)
        : rmpcc_transport_reference_tangent(e0, currentTau)
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
    _diagnostics.commandedLinearTwistNorm =
        _diagnostics.bodyTwist.head<3>().norm();
    _diagnostics.commandedAngularTwistNorm =
        _diagnostics.bodyTwist.tail<3>().norm();
    const Eigen::Vector<double,TWIST_DIM> pathVelocityResidual =
        _diagnostics.bodyTwist - _diagnostics.feedforwardBodyTwist;
    _diagnostics.pathVelocityLinearResidual =
        pathVelocityResidual.head<3>().norm();
    _diagnostics.pathVelocityAngularResidual =
        pathVelocityResidual.tail<3>().norm();
    const Eigen::Vector<double,TWIST_DIM> currentLagError = associatedPhaseFamily
        ? currentAssociatedResidual.vectorLag
        : (phaseFactorialFamily ? currentPhaseResidual.vectorLag
                                : currentLagProjection * currentCostError);
    const Eigen::Matrix<double, 1, TWIST_DIM> lagRow =
        (currentErrorTangent.transpose() * _rmpcc.metric) / std::sqrt(currentDenom);
    _diagnostics.lagError = (lagRow * e0)(0);
    _diagnostics.lagErrorNorm = currentLagError.norm();
    _diagnostics.lagTranslationErrorNorm = currentLagError.head<3>().norm();
    _diagnostics.lagRotationErrorNorm = currentLagError.tail<3>().norm();
    _diagnostics.contourError = currentContourError;
    _diagnostics.metricPhaseCorrection = currentPhaseResidual.metricPhaseCorrection;
    _diagnostics.taskPhaseCorrection = currentPhaseResidual.taskPhaseCorrection;
    _diagnostics.phaseContamination = currentPhaseResidual.metricPhaseCorrection
        - currentPhaseResidual.taskPhaseCorrection;
    _diagnostics.activePhaseCorrection = associatedPhaseFamily
        ? currentAssociatedResidual.context.phaseCorrection
        : currentPhaseResidual.phaseCorrection;
    _diagnostics.scalarLag = associatedPhaseFamily
        ? currentAssociatedResidual.context.positionLag
        : currentPhaseResidual.scalarLag;
    _diagnostics.poseArcLag = associatedPhaseFamily
        ? currentAssociatedResidual.scalarPoseArcLag : 0.0;
    _diagnostics.associatedProgress = associatedPhaseFamily
        ? currentAssociatedResidual.context.associatedProgress : _pathProgress;
    _diagnostics.phaseDenominator = associatedPhaseFamily
        ? currentAssociatedResidual.context.phaseDenominator
        : currentPhaseResidual.phaseDenominator;
    _diagnostics.phaseObservable = associatedPhaseFamily
        ? currentAssociatedResidual.context.observable
        : currentPhaseResidual.phaseObservable;
    _diagnostics.contourResidualNorm = associatedPhaseFamily
        ? currentAssociatedResidual.contour.norm()
        : currentPhaseResidual.contour.norm();
    _diagnostics.associatedContourResidualNorm = associatedPhaseFamily
        ? currentAssociatedResidual.contour.norm() : 0.0;
    _diagnostics.poseArcTableResolution =
        static_cast<double>(_poseArcTable->sample_count());
    _diagnostics.poseArcTableMinimumDensity =
        _poseArcTable->sample_count() > 0 ? _poseArcTable->minimum_density() : 0.0;
    _diagnostics.poseArcTableMaximumDensity =
        _poseArcTable->sample_count() > 0 ? _poseArcTable->maximum_density() : 0.0;
    _diagnostics.poseArcTableTotalLength =
        _poseArcTable->sample_count() > 0 ? _poseArcTable->total_arc_length() : 0.0;
    _diagnostics.trackingCostScale = _rmpcc.trackingCostScale;
    _diagnostics.vectorLagTranslationNorm =
        currentPhaseResidual.vectorLag.head<3>().norm();
    _diagnostics.vectorLagRotationNorm =
        currentPhaseResidual.vectorLag.tail<3>().norm();
    _diagnostics.parentFrameMotionActive = parentMotionActive;
    _diagnostics.parentFrameBodyTwist = _parentFrameMotion.body_twist();
    _diagnostics.measuredParentPose = _parentFrameMotion.current_pose();
    _diagnostics.predictedParentPoseFirst = parentTransform(1);
    _diagnostics.predictedParentPoseHorizon = parentTransform(N);
    const Eigen::Matrix4d firstPathPose = referenceTransform(_pathProgress);
    _diagnostics.parentReferenceFactorFirst = parent_frame_reference_factor(
        firstPathPose, parentTransform(0), parentTransform(1));
    _diagnostics.repairedReferenceDisplacementFirst =
        legacy_repaired_reference_displacement(
            firstPathPose, currentTau, _diagnostics.progressRate, dt,
            parentTransform(0), parentTransform(1));

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
    const VectorXd predictedScalarLagResiduals =
        scalarLagResidualSensitivity * zOpt + scalarLagResidualOffset;
    const VectorXd predictedPathVelocityResiduals =
        pathVelocityResidualSensitivity * zOpt + pathVelocityResidualOffset;
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
        const double lagCost = usesScalarLag[static_cast<size_t>(stage)]
            ? scalarLagCostWeights[static_cast<size_t>(stage)]
                * predictedScalarLagResiduals(stage)
                * predictedScalarLagResiduals(stage)
            : (lagError.transpose()
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

        Eigen::Vector<double,TWIST_DIM> residual;
        if(_rmpcc.referenceMotion == RmpccReferenceMotion::StageConsistent
           || parentRepairActive)
        {
            residual = predictedPathVelocityResiduals.segment<TWIST_DIM>(row);
        }
        else
        {
            residual = zOpt.segment<TWIST_DIM>(uOffset)
                - transportedTangents[static_cast<size_t>(stage)] * zOpt(sOffset);
        }
        _diagnostics.runningPathVelocityCost +=
            (residual.transpose() * pathVelocityWeight * residual)(0);
    }
    _diagnostics.referencePoseQueryTimeSeconds = referencePoseQueryTime;
    _diagnostics.referenceTangentQueryTimeSeconds = referenceTangentQueryTime;
    _diagnostics.referencePoseQueryCount = referencePoseQueryCount;
    _diagnostics.referenceTangentQueryCount = referenceTangentQueryCount;
    _diagnostics.referencePoseRequestCount = referencePoseRequestCount;
    _diagnostics.referenceTangentRequestCount = referenceTangentRequestCount;
    _diagnostics.stageRolloutLinearizationTimeSeconds = rolloutLinearizationTime;
    _diagnostics.residualHessianAssemblyTimeSeconds = residualHessianTime;
    _diagnostics.totalSolveTimeSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - totalSolveStart).count();
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
    _diagnostics.realizedLinearTwistNorm = realizedBodyTwist.head<3>().norm();
    _diagnostics.realizedAngularTwistNorm = realizedBodyTwist.tail<3>().norm();
    _diagnostics.linearTwistRealizationError =
        (realizedBodyTwist.head<3>() - _diagnostics.bodyTwist.head<3>()).norm();
    _diagnostics.angularTwistRealizationError =
        (realizedBodyTwist.tail<3>() - _diagnostics.bodyTwist.tail<3>()).norm();
    _diagnostics.externalLinearSaturationActive =
        _diagnostics.linearTwistRealizationError > 1e-6 ? 1.0 : 0.0;
    _diagnostics.externalAngularSaturationActive =
        _diagnostics.angularTwistRealizationError > 1e-6 ? 1.0 : 0.0;

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
    const bool parentMotionActive = _parentFrameMotion.has_velocity()
        && _parentFrameMotion.body_twist().squaredNorm() > 0.0;
    if(parentMotionActive
       && _rmpcc.referenceMotion == RmpccReferenceMotion::LegacyTangentProduct)
    {
        const auto referenceTangent = [this](const double progress)
        {
            return _trajectory.tangent_at_progress(progress, _rmpcc.tangentStep);
        };
        const auto parentTransform = [this, dt](const int stage)
        {
            return _parentFrameMotion.predicted_pose(stage, dt);
        };
        _predictedNextError = rmpcc_parent_repaired_legacy_state_step(
            realizedState, realizedInput, dt, 0,
            referenceTransform, referenceTangent, parentTransform).head<TWIST_DIM>();
    }
    else if(_rmpcc.predictorGeometry == RmpccPredictorGeometry::ExactSE3)
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
