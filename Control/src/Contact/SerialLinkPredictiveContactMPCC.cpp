/** @file SerialLinkPredictiveContactMPCC.cpp */

#include <Control/Contact/SerialLinkPredictiveContactMPCC.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace RobotLibrary { namespace Control {

void
SerialLinkPredictiveContactMPCC::set_predictive_contact_parameters(
    const PredictiveContactMpccParameters &parameters)
{
    const double values[] = {
        parameters.identifiedDt, parameters.stiffness, parameters.targetForce,
        parameters.forceWeight, parameters.minimumForce, parameters.maximumForce,
        parameters.maximumPenetrationIncrement, parameters.forceSlackWeight,
        parameters.penetrationSlackWeight, parameters.maximumForceSlack,
        parameters.maximumPenetrationSlack, parameters.tangentPositionWeight,
        parameters.pathLagPositionWeight, parameters.realizationAutoregressive,
        parameters.realizationInputGain, parameters.minimumContactModelForce,
        parameters.maximumRobotNormalCommand,
        parameters.maximumRobotNormalCommandStep,
        parameters.normalCommandSmoothWeight};
    for(const double value : values)
    {
        if(not std::isfinite(value) || value < 0.0)
        {
            throw std::invalid_argument(
                "Predictive contact MPCC parameters must be finite and nonnegative.");
        }
    }
    if(parameters.identifiedDt <= 0.0 || parameters.stiffness <= 0.0
       || parameters.targetForce <= 0.0
       || parameters.maximumPenetrationIncrement <= 0.0
       || parameters.minimumForce > parameters.maximumForce
       || parameters.realizationAutoregressive >= 1.0
       || parameters.realizationDelay < 0
       || parameters.maximumRobotNormalCommand <= 0.0
       || parameters.maximumRobotNormalCommandStep <= 0.0
       || not parameters.compressionDirectionParent.allFinite()
       || parameters.compressionDirectionParent.norm() <= 1e-12
       || not parameters.contactOffsetEndpoint.allFinite())
    {
        throw std::invalid_argument("Invalid predictive contact MPCC parameters.");
    }
    _parameters = parameters;
    _parameters.compressionDirectionParent.normalize();
    _diagnostics = PredictiveContactMpccDiagnostics{};
    _pastRobotNormalCommands.assign(
        static_cast<size_t>(_parameters.realizationDelay), 0.0);
}

void
SerialLinkPredictiveContactMPCC::set_force_measurement(
    const double measuredForce, const bool valid)
{
    if(valid && (not std::isfinite(measuredForce) || measuredForce < 0.0))
    {
        throw std::invalid_argument(
            "A valid predictive contact-force measurement must be finite and nonnegative.");
    }
    _measurementValid = valid;
    _measuredForce = valid ? measuredForce : 0.0;
    _contactModelValid = valid
        && measuredForce >= _parameters.minimumContactModelForce;
}

void
SerialLinkPredictiveContactMPCC::set_normal_realization_state(
    const double measuredRobotNormalVelocity,
    const std::vector<double> &pastRobotNormalCommands)
{
    if(not std::isfinite(measuredRobotNormalVelocity)
       || static_cast<int>(pastRobotNormalCommands.size())
            != _parameters.realizationDelay)
    {
        throw std::invalid_argument("Invalid measured realization state/history.");
    }
    for(const double value : pastRobotNormalCommands)
    {
        if(not std::isfinite(value))
            throw std::invalid_argument("Issued normal-command history must be finite.");
    }
    _measuredRobotNormalVelocity = measuredRobotNormalVelocity;
    _pastRobotNormalCommands = pastRobotNormalCommands;
}

bool
SerialLinkPredictiveContactMPCC::prediction_enabled() const
{
    return _parameters.mode != PredictiveContactMode::Disabled
        && _contactModelValid;
}

PredictiveContactKinematics
SerialLinkPredictiveContactMPCC::make_kinematics(
    const MpccQpExtensionContext &context) const
{
    PredictiveContactKinematics kinematics;
    kinematics.horizon = context.horizon;
    kinematics.stageControlDimension = context.stageControlDimension;
    kinematics.dt = context.dt;
    kinematics.expectedDt = _parameters.identifiedDt;
    kinematics.predictionRotation = context.predictionRotation;
    kinematics.endpointRotationBase = context.currentToolRotationBase;
    kinematics.endpointPositionBase = context.currentToolPositionBase;
    kinematics.contactOffsetEndpoint = _parameters.contactOffsetEndpoint;
    kinematics.compressionDirectionParent =
        _parameters.compressionDirectionParent;
    kinematics.parentTransforms = context.parentTransforms;
    kinematics.actuationAware = _parameters.actuationAware;
    kinematics.realizationAutoregressive =
        _parameters.realizationAutoregressive;
    kinematics.realizationInputGain = _parameters.realizationInputGain;
    kinematics.realizationDelay = _parameters.realizationDelay;
    kinematics.initialRealizedRobotNormalVelocity =
        _measuredRobotNormalVelocity;
    kinematics.pastRobotNormalCommands = _pastRobotNormalCommands;
    return kinematics;
}

Eigen::Matrix3d
SerialLinkPredictiveContactMPCC::position_error_weight(
    const int stage,
    const Eigen::Vector<double,ERROR_DIM> &pathTangent,
    const Eigen::Matrix3d &defaultWeight,
    const Eigen::Matrix3d &predictionRotation,
    const Eigen::Matrix4d &predictedParentTransform) const
{
    (void)stage;
    if(_parameters.mode != PredictiveContactMode::Active
       || not prediction_enabled())
    {
        return defaultWeight;
    }
    const Eigen::Vector3d normal = (
        predictionRotation.transpose()
        * predictedParentTransform.block<3,3>(0,0)
        * _parameters.compressionDirectionParent).normalized();
    Eigen::Vector3d tangent = pathTangent.head<3>();
    tangent -= normal * normal.dot(tangent);
    if(tangent.norm() <= 1e-10)
    {
        tangent = normal.unitOrthogonal();
    }
    tangent.normalize();
    const Eigen::Vector3d crossTrack = normal.cross(tangent).normalized();
    return _parameters.pathLagPositionWeight * tangent * tangent.transpose()
        + _parameters.tangentPositionWeight
            * crossTrack * crossTrack.transpose();
}

void
SerialLinkPredictiveContactMPCC::extend_qp_problem(
    const MpccQpExtensionContext &context,
    Eigen::MatrixXd &hessian,
    Eigen::VectorXd &gradient,
    Eigen::MatrixXd &constraintMatrix,
    Eigen::VectorXd &constraintVector,
    Eigen::VectorXd &seed)
{
    _diagnostics = PredictiveContactMpccDiagnostics{};
    _diagnostics.mode = _parameters.mode;
    _diagnostics.forceValid = _measurementValid;
    _diagnostics.contactModelValid = _contactModelValid;
    _diagnostics.measuredForce = _measuredForce;
    _diagnostics.stiffness = _parameters.stiffness;
    _diagnostics.measuredRobotNormalVelocity =
        _measuredRobotNormalVelocity;
    _baseControlDimension = context.baseControlDimension;
    if(not prediction_enabled())
    {
        return;
    }
    if(context.stageControlDimension != 7
       || context.baseControlDimension != 7 * context.horizon)
    {
        throw std::runtime_error("CA-MPCC-C requires exactly seven inputs per stage.");
    }

    _lastModel = build_predictive_contact_affine_model(make_kinematics(context));
    _lastForceMap = _parameters.stiffness
        * _lastModel.penetrationIncrementMap;
    _lastForceConstant = Eigen::VectorXd::Constant(
        context.horizon, _measuredForce)
        + _parameters.stiffness * _lastModel.penetrationIncrementOffset;
    _diagnostics.forceMap = _lastForceMap;
    _diagnostics.forceOffset = _parameters.stiffness
        * _lastModel.penetrationIncrementOffset;

    if(_parameters.mode == PredictiveContactMode::Shadow)
    {
        return;
    }

    const int N = context.horizon;
    const int C = context.baseControlDimension;
    const int V = C + 2 * N;
    Eigen::MatrixXd augmentedHessian = Eigen::MatrixXd::Zero(V, V);
    augmentedHessian.topLeftCorner(C, C) = hessian;
    hessian = std::move(augmentedHessian);
    Eigen::VectorXd augmentedGradient = Eigen::VectorXd::Zero(V);
    augmentedGradient.head(C) = gradient;
    gradient = std::move(augmentedGradient);
    Eigen::MatrixXd augmentedConstraints = Eigen::MatrixXd::Zero(
        constraintMatrix.rows(), V);
    augmentedConstraints.leftCols(C) = constraintMatrix;
    constraintMatrix = std::move(augmentedConstraints);
    Eigen::VectorXd augmentedSeed = Eigen::VectorXd::Zero(V);
    augmentedSeed.head(C) = seed;
    if(context.previousWarmStart.size() == V)
    {
        augmentedSeed.tail(2 * N) =
            context.previousWarmStart.tail(2 * N).cwiseMax(0.0);
    }
    seed = std::move(augmentedSeed);

    Eigen::MatrixXd forceMap = Eigen::MatrixXd::Zero(N, V);
    forceMap.leftCols(C) = _lastForceMap;
    const Eigen::VectorXd forceGradientBefore = gradient;
    add_normalized_predictive_force_cost(
        forceMap, _lastForceConstant, _parameters.targetForce,
        _parameters.forceWeight, hessian, gradient);
    _diagnostics.forceGradientNorm =
        (gradient - forceGradientBefore).head(C).norm();
    _diagnostics.forceObjectiveActive = _parameters.forceWeight > 0.0;

    // Keep the commanded robot-side normal action in the E06-Q force-safe
    // small-signal domain and suppress horizon-to-horizon high-frequency action.
    const Eigen::MatrixXd &commandMap =
        _lastModel.commandedRobotNormalVelocityMap;
    Eigen::MatrixXd slewMap = Eigen::MatrixXd::Zero(N, C);
    Eigen::VectorXd slewOffset = Eigen::VectorXd::Zero(N);
    slewMap.row(0) = commandMap.row(0);
    const double previousCommand = _pastRobotNormalCommands.empty()
        ? 0.0 : _pastRobotNormalCommands.back();
    slewOffset(0) = -previousCommand;
    for(int stage = 1; stage < N; ++stage)
        slewMap.row(stage) = commandMap.row(stage) - commandMap.row(stage - 1);
    if(_parameters.normalActionGuardEnabled)
    {
        const double smoothScale = _parameters.normalCommandSmoothWeight
            / (_parameters.maximumRobotNormalCommand
               * _parameters.maximumRobotNormalCommand);
        hessian.topLeftCorner(C, C) +=
            2.0 * smoothScale * slewMap.transpose() * slewMap;
        gradient.head(C) +=
            2.0 * smoothScale * slewMap.transpose() * slewOffset;

        // The base MPCC warm start is feasible for its own box/progress
        // constraints, but it predates the contact-normal magnitude and slew
        // rows below.  Active-set QP requires a feasible seed, so project each
        // stage's scalar normal action into the frozen guard before adding the
        // soft-constraint seeds.  This changes only the numerical seed, never
        // the objective or feasible set.
        const double heldSeedCommand = std::clamp(
            previousCommand,
            -_parameters.maximumRobotNormalCommand,
            _parameters.maximumRobotNormalCommand);
        for(int stage = 0; stage < N; ++stage)
        {
            const Eigen::RowVectorXd row = commandMap.row(stage);
            const double requested = row.dot(seed.head(C));
            const double normSquared = row.squaredNorm();
            if(normSquared > 1e-18)
            {
                seed.head(C) += row.transpose()
                    * ((heldSeedCommand - requested) / normSquared);
            }
        }
        if((constraintMatrix * seed - constraintVector).maxCoeff() > 1e-9)
        {
            for(int stage = 0; stage < N; ++stage)
                seed.segment(stage * context.stageControlDimension, 6).setZero();
            for(int stage = 0; stage < N; ++stage)
            {
                const Eigen::RowVectorXd row = commandMap.row(stage);
                const double normSquared = row.squaredNorm();
                if(normSquared > 1e-18)
                    seed.head(C) += row.transpose()
                        * (heldSeedCommand / normSquared);
            }
        }
    }
    for(int stage = 0; stage < N; ++stage)
    {
        hessian(C + stage, C + stage) +=
            2.0 * _parameters.forceSlackWeight;
        hessian(C + N + stage, C + N + stage) +=
            2.0 * _parameters.penetrationSlackWeight
                * _parameters.maximumPenetrationIncrement
                * _parameters.maximumPenetrationIncrement;
    }

    const int oldRows = constraintMatrix.rows();
    const int actionRows = _parameters.normalActionGuardEnabled ? 4 * N : 0;
    // Force and penetration bands are genuinely soft: each slack is
    // constrained nonnegative but deliberately has no hard upper bound.
    // A capped slack can make a model-validity guard infeasible precisely
    // when its diagnostic violation is largest.
    Eigen::MatrixXd bounded = Eigen::MatrixXd::Zero(
        oldRows + 6 * N + actionRows, V);
    bounded.topRows(oldRows) = constraintMatrix;
    Eigen::VectorXd limits = Eigen::VectorXd::Zero(
        oldRows + 6 * N + actionRows);
    limits.head(oldRows) = constraintVector;

    for(int stage = 0; stage < N; ++stage)
    {
        const int row = oldRows + 6 * stage;
        const int actionRow = oldRows + 6 * N + 4 * stage;
        const int forceSlack = C + stage;
        const int penetrationSlack = C + N + stage;
        const Eigen::RowVectorXd forceRow = _lastForceMap.row(stage);
        const double forceConstant = _lastForceConstant(stage);
        const Eigen::RowVectorXd penetrationRow =
            _lastModel.penetrationIncrementMap.row(stage);
        const double penetrationConstant =
            _lastModel.penetrationIncrementOffset(stage);

        bounded.block(row, 0, 1, C) = forceRow;
        bounded(row, forceSlack) = -1.0;
        limits(row) = _parameters.maximumForce - forceConstant;
        bounded.block(row + 1, 0, 1, C) = -forceRow;
        bounded(row + 1, forceSlack) = -1.0;
        limits(row + 1) = forceConstant - _parameters.minimumForce;
        bounded(row + 2, forceSlack) = -1.0;
        limits(row + 2) = 0.0;
        // The penetration slack decision is dimensionless:
        // sigma_rho = epsilon_rho / Delta_rho_max.  Scaling this block avoids
        // putting metre-scale coefficients beside O(1) force constraints and
        // a 1e12--1e16 SI penalty in the same active-set KKT system.
        const double rhoScale = _parameters.maximumPenetrationIncrement;
        bounded.block(row + 3, 0, 1, C) = penetrationRow / rhoScale;
        bounded(row + 3, penetrationSlack) = -1.0;
        limits(row + 3) = 1.0 - penetrationConstant / rhoScale;
        bounded.block(row + 4, 0, 1, C) = -penetrationRow / rhoScale;
        bounded(row + 4, penetrationSlack) = -1.0;
        limits(row + 4) = 1.0 + penetrationConstant / rhoScale;
        bounded(row + 5, penetrationSlack) = -1.0;
        limits(row + 5) = 0.0;

        if(_parameters.normalActionGuardEnabled)
        {
            // Normalize both guards to O(0.1).  Their physical limits are
            // 1e-4 and 1e-6 m/s, respectively; leaving those raw beside the
            // force/penetration rows makes an otherwise identical feasible
            // set poorly resolved by the frozen outer-QP tolerance.
            const double commandScale =
                10.0 * _parameters.maximumRobotNormalCommand;
            const double slewScale =
                10.0 * _parameters.maximumRobotNormalCommandStep;
            bounded.block(actionRow, 0, 1, C) =
                commandMap.row(stage) / commandScale;
            limits(actionRow) =
                _parameters.maximumRobotNormalCommand / commandScale;
            bounded.block(actionRow + 1, 0, 1, C) =
                -commandMap.row(stage) / commandScale;
            limits(actionRow + 1) =
                _parameters.maximumRobotNormalCommand / commandScale;
            bounded.block(actionRow + 2, 0, 1, C) =
                slewMap.row(stage) / slewScale;
            limits(actionRow + 2) =
                _parameters.maximumRobotNormalCommandStep / slewScale
                - slewOffset(stage) / slewScale;
            bounded.block(actionRow + 3, 0, 1, C) =
                -slewMap.row(stage) / slewScale;
            limits(actionRow + 3) =
                _parameters.maximumRobotNormalCommandStep / slewScale
                + slewOffset(stage) / slewScale;
        }

        const double seededForce = forceConstant
            + forceRow.dot(seed.head(C));
        seed(forceSlack) = std::max({
            seed(forceSlack), 0.0,
            seededForce - _parameters.maximumForce,
            _parameters.minimumForce - seededForce});
        const double seededPenetration = penetrationConstant
            + penetrationRow.dot(seed.head(C));
        seed(penetrationSlack) = std::max(
            seed(penetrationSlack),
            std::max(0.0, std::abs(seededPenetration) / rhoScale - 1.0));
    }
    constraintMatrix = std::move(bounded);
    constraintVector = std::move(limits);
}

void
SerialLinkPredictiveContactMPCC::shift_extension_warm_start(
    const MpccQpExtensionContext &context,
    const Eigen::VectorXd &optimum,
    Eigen::VectorXd &shiftedWarmStart)
{
    if(_parameters.mode != PredictiveContactMode::Active
       || not prediction_enabled())
    {
        return;
    }
    const int C = context.baseControlDimension;
    const int N = context.horizon;
    if(optimum.size() != C + 2 * N || shiftedWarmStart.size() != C + 2 * N)
    {
        return;
    }
    for(int stage = 0; stage < N - 1; ++stage)
    {
        shiftedWarmStart(C + stage) = optimum(C + stage + 1);
        shiftedWarmStart(C + N + stage) = optimum(C + N + stage + 1);
    }
}

void
SerialLinkPredictiveContactMPCC::on_extended_qp_solution(
    const MpccQpExtensionContext &context,
    const Eigen::VectorXd &optimum)
{
    if(not prediction_enabled())
    {
        return;
    }
    const Eigen::VectorXd decision = optimum.head(context.baseControlDimension);
    _diagnostics.predictedRelativeNormalVelocity =
        _lastModel.relativeVelocityMap * decision
        + _lastModel.relativeVelocityOffset;
    _diagnostics.predictedPenetrationIncrement =
        _lastModel.penetrationIncrementMap * decision
        + _lastModel.penetrationIncrementOffset;
    _diagnostics.predictedForce =
        _lastForceConstant + _lastForceMap * decision;
    _diagnostics.predictedCommandedRobotNormalVelocity =
        _lastModel.commandedRobotNormalVelocityMap * decision;
    _diagnostics.predictedRealizedRobotNormalVelocity =
        _lastModel.realizedRobotNormalVelocityMap * decision
        + _lastModel.realizedRobotNormalVelocityOffset;
    _diagnostics.optimizedRelativeNormalVelocityStage0 =
        _diagnostics.predictedRelativeNormalVelocity(0);
    _diagnostics.optimizedRobotNormalCommandStage0 =
        _diagnostics.predictedCommandedRobotNormalVelocity(0);
    const double previousCommand = _pastRobotNormalCommands.empty()
        ? 0.0 : _pastRobotNormalCommands.back();
    _diagnostics.robotNormalCommandSlewStage0 =
        _diagnostics.optimizedRobotNormalCommandStage0 - previousCommand;
    Eigen::VectorXd commandSlew =
        _diagnostics.predictedCommandedRobotNormalVelocity;
    if(commandSlew.size() > 0)
    {
        for(int stage = commandSlew.size() - 1; stage > 0; --stage)
            commandSlew(stage) -= commandSlew(stage - 1);
        commandSlew(0) -= previousCommand;
        _diagnostics.normalCommandSmoothCost =
            _parameters.normalCommandSmoothWeight
            * commandSlew.squaredNorm()
            / (_parameters.maximumRobotNormalCommand
               * _parameters.maximumRobotNormalCommand);
    }
    const Eigen::VectorXd normalizedResidual =
        (_diagnostics.predictedForce.array() - _parameters.targetForce)
        / _parameters.targetForce;
    _diagnostics.forceCost =
        _parameters.forceWeight * normalizedResidual.squaredNorm();
    if(_parameters.mode == PredictiveContactMode::Active
       && optimum.size() == context.baseControlDimension + 2 * context.horizon)
    {
        _diagnostics.maximumForceSlack = std::max(
            0.0, optimum.segment(
                context.baseControlDimension, context.horizon).maxCoeff());
        _diagnostics.maximumPenetrationSlack = std::max(
            0.0, _parameters.maximumPenetrationIncrement
                * optimum.tail(context.horizon).maxCoeff());
    }
}

void
SerialLinkPredictiveContactMPCC::reset_additional_virtual_state()
{
    _diagnostics = PredictiveContactMpccDiagnostics{};
    _lastForceMap.resize(0, 0);
    _lastForceConstant.resize(0);
    _baseControlDimension = 0;
    _measurementValid = false;
    _contactModelValid = false;
    _measuredRobotNormalVelocity = 0.0;
    _pastRobotNormalCommands.assign(
        static_cast<size_t>(_parameters.realizationDelay), 0.0);
}

} } // namespace RobotLibrary::Control
