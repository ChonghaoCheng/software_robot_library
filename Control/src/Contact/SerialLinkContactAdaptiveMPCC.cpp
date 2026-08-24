/** @file SerialLinkContactAdaptiveMPCC.cpp */

#include <Control/Contact/SerialLinkContactAdaptiveMPCC.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace RobotLibrary { namespace Control {

double
contact_force_deadband(const double forceError, const double deadband)
{
    if(not std::isfinite(forceError) or not std::isfinite(deadband) or deadband < 0.0)
    {
        throw std::invalid_argument("Contact-force deadband inputs must be finite and nonnegative.");
    }
    if(forceError > deadband) return forceError - deadband;
    if(forceError < -deadband) return forceError + deadband;
    return 0.0;
}

double
contact_preferred_delta_rate(
    const ContactAdaptiveMpccParameters &parameters,
    const double measuredForce)
{
    if(not std::isfinite(measuredForce))
    {
        throw std::invalid_argument("Measured normal force must be finite.");
    }
    if(not parameters.adaptationEnabled)
    {
        return 0.0;
    }
    const double forceError = parameters.targetForce - measuredForce;
    return std::clamp(
        parameters.deltaGain
            * contact_force_deadband(forceError, parameters.forceDeadband),
        -parameters.deltaRateMaximum,
        parameters.deltaRateMaximum);
}

Eigen::Vector<double,6>
contact_delta_tangent_body(
    const Eigen::Matrix3d &pathRotationInParent,
    const Eigen::Vector3d &compressionDirectionInParent)
{
    Eigen::Vector<double,6> tangent = Eigen::Vector<double,6>::Zero();
    tangent.head<3>() = pathRotationInParent.transpose()
        * compressionDirectionInParent.normalized();
    return tangent;
}

Eigen::Vector<double,6>
contact_delta_tangent_prediction_frame(
    const Eigen::Matrix3d &predictionRotationInBase,
    const Eigen::Matrix3d &parentRotationInBase,
    const Eigen::Vector3d &compressionDirectionInParent)
{
    Eigen::Vector<double,6> tangent = Eigen::Vector<double,6>::Zero();
    tangent.head<3>() = predictionRotationInBase.transpose()
        * parentRotationInBase * compressionDirectionInParent.normalized();
    return tangent;
}

void
SerialLinkContactAdaptiveMPCC::set_contact_adaptive_parameters(
    const ContactAdaptiveMpccParameters &parameters)
{
    if(not std::isfinite(parameters.targetForce) or parameters.targetForce < 0.0
       or not std::isfinite(parameters.deltaGain) or parameters.deltaGain < 0.0
       or not std::isfinite(parameters.forceDeadband) or parameters.forceDeadband < 0.0
       or not std::isfinite(parameters.deltaMinimum)
       or not std::isfinite(parameters.deltaMaximum)
       or parameters.deltaMinimum > parameters.deltaMaximum
       or not std::isfinite(parameters.deltaRateMaximum)
       or parameters.deltaRateMaximum < 0.0
       or not std::isfinite(parameters.forceRateWeight)
       or parameters.forceRateWeight < 0.0
       or not std::isfinite(parameters.deltaRateSmoothWeight)
       or parameters.deltaRateSmoothWeight < 0.0
       or not std::isfinite(parameters.normalAdmittanceGain)
       or parameters.normalAdmittanceGain < 0.0
       or not parameters.compressionDirectionInParent.allFinite()
       or parameters.compressionDirectionInParent.norm() <= 1e-12)
    {
        throw std::invalid_argument("Invalid contact-adaptive MPCC parameters.");
    }
    _parameters = parameters;
    _parameters.compressionDirectionInParent.normalize();
    _delta = std::clamp(_delta, _parameters.deltaMinimum, _parameters.deltaMaximum);
}

void
SerialLinkContactAdaptiveMPCC::set_measured_normal_force(const double measuredForce)
{
    if(not std::isfinite(measuredForce))
    {
        throw std::invalid_argument("Measured normal force must be finite.");
    }
    _measuredForce = measuredForce;
}

void
SerialLinkContactAdaptiveMPCC::set_contact_coordinate(const double delta)
{
    if(not std::isfinite(delta)
       or delta < _parameters.deltaMinimum
       or delta > _parameters.deltaMaximum)
    {
        throw std::invalid_argument("Contact coordinate is outside configured bounds.");
    }
    _delta = delta;
}

void
SerialLinkContactAdaptiveMPCC::update_contact_coordinate_from_force(const double dt)
{
    if(not std::isfinite(dt) or dt <= 0.0)
    {
        throw std::invalid_argument("Contact-coordinate update dt must be positive.");
    }
    if(not _parameters.adaptationEnabled or _parameters.optimizeDeltaRate)
    {
        _lastDeltaRate = 0.0;
        return;
    }
    const double preferred = contact_preferred_delta_rate(_parameters, _measuredForce);
    const double nextDelta = std::clamp(
        _delta + dt * preferred,
        _parameters.deltaMinimum, _parameters.deltaMaximum);
    _lastDeltaRate = (nextDelta - _delta) / dt;
    _delta = nextDelta;
}

RobotLibrary::Model::Pose
SerialLinkContactAdaptiveMPCC::reference_pose_at_progress(const double progress)
{
    RobotLibrary::Model::Pose poseInParent =
        reference_trajectory().pose_at_progress(progress);
    if(_parameters.adaptationEnabled)
    {
        poseInParent = RobotLibrary::Model::Pose(
            poseInParent.translation()
                + _delta * _parameters.compressionDirectionInParent,
            poseInParent.quaternion());
    }
    return RobotLibrary::Trajectory::express_pose_in_base(
        trajectory_frame(), poseInParent);
}

Eigen::MatrixXd
SerialLinkContactAdaptiveMPCC::additional_reference_tangents(
    const int stage,
    const double progress,
    const Eigen::Matrix3d &predictionRotation,
    const Eigen::Matrix4d &predictedParentTransform) const
{
    (void)stage;
    (void)progress;
    const int additionalCount =
        stage_control_dimension() - BASE_STAGE_CONTROL_DIMENSION;
    Eigen::MatrixXd tangent = Eigen::MatrixXd::Zero(6, additionalCount);
    if(_parameters.adaptationEnabled and _parameters.optimizeDeltaRate)
    {
        tangent.col(0) = contact_delta_tangent_prediction_frame(
            predictionRotation,
            predictedParentTransform.block<3,3>(0,0),
            _parameters.compressionDirectionInParent);
    }
    return tangent;
}

void
SerialLinkContactAdaptiveMPCC::configure_additional_stage_inputs(
    const int stage,
    Eigen::VectorXd &lower,
    Eigen::VectorXd &upper,
    Eigen::VectorXd &nominal) const
{
    const int index = stage * stage_control_dimension() + DELTA_RATE_INDEX;
    if(not _parameters.optimizeDeltaRate)
    {
        return;
    }
    const double limit = _parameters.adaptationEnabled
        ? _parameters.deltaRateMaximum : 0.0;
    lower(index) = -limit;
    upper(index) = limit;
    nominal(index) = contact_preferred_delta_rate(_parameters, _measuredForce);
}

void
SerialLinkContactAdaptiveMPCC::extend_qp_problem(
    const MpccQpExtensionContext &context,
    Eigen::MatrixXd &hessian,
    Eigen::VectorXd &gradient,
    Eigen::MatrixXd &constraintMatrix,
    Eigen::VectorXd &constraintVector,
    Eigen::VectorXd &seed)
{
    const int N = context.horizon;
    const int D = context.stageControlDimension;
    if(D == 7 and not _parameters.optimizeDeltaRate)
    {
        return;
    }
    if(D != 8 or context.baseControlDimension != N * D)
    {
        throw std::runtime_error("CA-MPCC-B requires an interleaved 8D stage decision.");
    }

    const double preferred = contact_preferred_delta_rate(_parameters, _measuredForce);
    if(_parameters.adaptationEnabled and _parameters.contactCostEnabled)
    {
        for(int stage = 0; stage < N; ++stage)
        {
            const int index = stage * D + DELTA_RATE_INDEX;
            hessian(index, index) += _parameters.forceRateWeight;
            gradient(index) -= _parameters.forceRateWeight * preferred;
        }

        if(_parameters.deltaRateSmoothWeight > 0.0)
        {
            Eigen::MatrixXd difference = Eigen::MatrixXd::Zero(N, N * D);
            Eigen::VectorXd previous = Eigen::VectorXd::Zero(N);
            for(int stage = 0; stage < N; ++stage)
            {
                difference(stage, stage * D + DELTA_RATE_INDEX) = 1.0;
                if(stage > 0)
                {
                    difference(stage, (stage - 1) * D + DELTA_RATE_INDEX) = -1.0;
                }
            }
            previous(0) = _lastDeltaRate;
            hessian += _parameters.deltaRateSmoothWeight
                * difference.transpose() * difference;
            gradient -= _parameters.deltaRateSmoothWeight
                * difference.transpose() * previous;
        }
    }

    const int oldRows = constraintMatrix.rows();
    Eigen::MatrixXd bounded = Eigen::MatrixXd::Zero(oldRows + 2 * N, N * D);
    bounded.topRows(oldRows) = constraintMatrix;
    Eigen::VectorXd limits = Eigen::VectorXd::Zero(oldRows + 2 * N);
    limits.head(oldRows) = constraintVector;
    for(int state = 1; state <= N; ++state)
    {
        const int upperRow = oldRows + 2 * (state - 1);
        const int lowerRow = upperRow + 1;
        for(int input = 0; input < state; ++input)
        {
            const int index = input * D + DELTA_RATE_INDEX;
            bounded(upperRow, index) = context.dt;
            bounded(lowerRow, index) = -context.dt;
        }
        limits(upperRow) = _parameters.deltaMaximum - _delta;
        limits(lowerRow) = _delta - _parameters.deltaMinimum;
    }
    constraintMatrix = std::move(bounded);
    constraintVector = std::move(limits);

    // Sequential projection preserves both the rate box and every cumulative
    // delta state, so the active-set solver always receives a feasible seed.
    double predictedDelta = _delta;
    for(int stage = 0; stage < N; ++stage)
    {
        const int index = stage * D + DELTA_RATE_INDEX;
        double rate = _parameters.adaptationEnabled
            ? std::clamp(seed(index), -_parameters.deltaRateMaximum,
                         _parameters.deltaRateMaximum)
            : 0.0;
        const double nextDelta = std::clamp(
            predictedDelta + context.dt * rate,
            _parameters.deltaMinimum, _parameters.deltaMaximum);
        seed(index) = (nextDelta - predictedDelta) / context.dt;
        predictedDelta = nextDelta;
    }
}

void
SerialLinkContactAdaptiveMPCC::on_extended_qp_solution(
    const MpccQpExtensionContext &context,
    const Eigen::VectorXd &optimum)
{
    const int N = context.horizon;
    const int D = context.stageControlDimension;
    const double preferred = contact_preferred_delta_rate(_parameters, _measuredForce);
    const double firstRate = _parameters.optimizeDeltaRate
        ? optimum(DELTA_RATE_INDEX) : _lastDeltaRate;
    const double nextDelta = std::clamp(
        _delta + context.dt * firstRate,
        _parameters.deltaMinimum, _parameters.deltaMaximum);

    _contactDiagnostics = ContactAdaptiveMpccDiagnostics{};
    _contactDiagnostics.measuredForce = _measuredForce;
    _contactDiagnostics.targetForce = _parameters.targetForce;
    _contactDiagnostics.forceError = _parameters.targetForce - _measuredForce;
    _contactDiagnostics.preferredDeltaRate = preferred;
    _contactDiagnostics.optimizedDeltaRate = firstRate;
    _contactDiagnostics.delta = _delta;
    _contactDiagnostics.nextDelta = nextDelta;
    _contactDiagnostics.deltaSaturated =
        std::abs(nextDelta - _parameters.deltaMinimum) <= 1e-10
        or std::abs(nextDelta - _parameters.deltaMaximum) <= 1e-10;
    _contactDiagnostics.deltaRateSaturated = _parameters.deltaRateMaximum > 0.0
        and std::abs(firstRate) >= _parameters.deltaRateMaximum - 1e-10;
    _contactDiagnostics.optimizedDeltaRates.resize(N);
    _contactDiagnostics.predictedDelta.resize(N);
    double predictedDelta = _delta;
    for(int stage = 0; stage < N; ++stage)
    {
        const double rate = _parameters.optimizeDeltaRate
            ? optimum(stage * D + DELTA_RATE_INDEX) : _lastDeltaRate;
        _contactDiagnostics.optimizedDeltaRates(stage) = rate;
        predictedDelta += context.dt * rate;
        _contactDiagnostics.predictedDelta(stage) = predictedDelta;
    }
    _contactDiagnostics.firstDeltaTangent =
        contact_delta_tangent_prediction_frame(
            context.predictionRotation,
            context.parentTransforms.front().block<3,3>(0,0),
            _parameters.compressionDirectionInParent);

    if(_parameters.optimizeDeltaRate)
    {
        _delta = nextDelta;
    }
    _lastDeltaRate = firstRate;
}

void
SerialLinkContactAdaptiveMPCC::reset_additional_virtual_state()
{
    _delta = 0.0;
    _lastDeltaRate = 0.0;
    _contactDiagnostics = ContactAdaptiveMpccDiagnostics{};
}

Eigen::Vector<double,6>
SerialLinkContactAdaptiveMPCC::postprocess_base_twist(
    const Eigen::Vector<double,6> &baseTwist,
    const double dt)
{
    (void)dt;
    if(not _parameters.normalAdmittanceEnabled)
    {
        return baseTwist;
    }
    const double normalRate = std::clamp(
        _parameters.normalAdmittanceGain
            * contact_force_deadband(
                _parameters.targetForce - _measuredForce,
                _parameters.forceDeadband),
        -_parameters.deltaRateMaximum,
        _parameters.deltaRateMaximum);
    const Eigen::Vector3d normalBase =
        trajectory_frame().transformInBase.block<3,3>(0,0)
        * _parameters.compressionDirectionInParent;
    Eigen::Vector<double,6> result = baseTwist;
    result.head<3>() += normalBase
        * (normalRate - normalBase.dot(result.head<3>()));
    return result;
}

} } // namespace RobotLibrary::Control
