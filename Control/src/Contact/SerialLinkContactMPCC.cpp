/** @file SerialLinkContactMPCC.cpp */

#include <Control/Contact/SerialLinkContactMPCC.h>

#include "detail/ContactGeometry.h"
#include "detail/ContactParameterValidation.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace RobotLibrary { namespace Control {

void
SerialLinkContactMPCC::set_contact_parameters(const ContactParameters &parameters)
{
    detail::validate_contact_parameters(parameters);
    _contactParameters = parameters;
    _contactDiagnostics = ContactMpcDiagnostics{};
}

void
SerialLinkContactMPCC::set_contact_state(const ContactState &state)
{
    if(!state.normalBase.allFinite() || !state.surfaceVelocityBase.allFinite())
    {
        throw std::invalid_argument("ContactState vectors must be finite.");
    }
    _contactState = state;
    _contactState.measuredNormalForce =
        detail::validate_measured_normal_force(state.measuredNormalForce);
}

bool
SerialLinkContactMPCC::contact_active() const
{
    return _contactParameters.mode != ContactMode::Disabled
        && (_contactState.inContact || _contactParameters.approachWhenNotInContact);
}

Eigen::Matrix3d
SerialLinkContactMPCC::position_error_weight(
    const int stage,
    const Eigen::Vector<double,ERROR_DIM> &pathTangent,
    const Eigen::Matrix3d &defaultWeight,
    const Eigen::Matrix3d &predictionRotation,
    const Eigen::Matrix4d &predictedParentTransform) const
{
    (void)stage;
    if(!contact_active())
    {
        return defaultWeight;
    }

    const Eigen::Matrix3d boardRotation =
        predictedParentTransform.block<3,3>(0,0);
    const Eigen::Vector3d normalBase = detail::inward_contact_normal(
        boardRotation, _contactParameters.normalAxisInBoard);
    const Eigen::Vector3d normal =
        (predictionRotation.transpose() * normalBase).normalized();

    Eigen::Vector3d tangent = pathTangent.head<3>();
    tangent -= normal * normal.dot(tangent);
    if(tangent.norm() <= 1e-9)
    {
        tangent = predictionRotation.transpose()
            * (boardRotation * _contactParameters.tangentAxisInBoard);
        tangent -= normal * normal.dot(tangent);
    }
    tangent.normalize();
    const Eigen::Vector3d crossTrack = normal.cross(tangent).normalized();

    return _contactParameters.pathLagPositionWeight * tangent * tangent.transpose()
        + _contactParameters.tangentPositionWeight * crossTrack * crossTrack.transpose()
        + _contactParameters.normalPositionWeight * normal * normal.transpose();
}

void
SerialLinkContactMPCC::extend_qp_problem(
    const MpccQpExtensionContext &context,
    Eigen::MatrixXd &hessian,
    Eigen::VectorXd &gradient,
    Eigen::MatrixXd &constraintMatrix,
    Eigen::VectorXd &constraintVector,
    Eigen::VectorXd &seed)
{
    using Eigen::MatrixXd;
    using Eigen::RowVectorXd;
    using Eigen::VectorXd;

    _contactDiagnostics = ContactMpcDiagnostics{};
    _contactDiagnostics.measuredNormalForce = _contactState.measuredNormalForce;
    _contactDiagnostics.contactModeActive = contact_active();
    if(!contact_active())
    {
        return;
    }

    const int N = context.horizon;
    const int C = context.baseControlDimension;
    const bool constraintMode = _contactParameters.mode == ContactMode::Constraint;
    const int V = C + (constraintMode ? N : 0);

    if(static_cast<int>(context.parentTransforms.size()) != N + 1)
    {
        throw std::runtime_error("Contact MPCC requires N+1 parent-frame predictions.");
    }

    const Eigen::Matrix4d &boardNow = context.parentTransforms.front();
    const Eigen::Vector3d normalNow = detail::inward_contact_normal(
        boardNow.block<3,3>(0,0), _contactParameters.normalAxisInBoard);
    _currentToolPosition = context.currentToolPositionBase;
    _currentSignedNormalCoordinate = detail::signed_normal_coordinate(
        normalNow, _currentToolPosition, boardNow.block<3,1>(0,3));
    _contactDiagnostics.currentSignedNormalCoordinate = _currentSignedNormalCoordinate;

    if(constraintMode)
    {
        MatrixXd augmentedH = MatrixXd::Zero(V, V);
        augmentedH.topLeftCorner(C, C) = hessian;
        hessian = std::move(augmentedH);

        VectorXd augmentedGradient = VectorXd::Zero(V);
        augmentedGradient.head(C) = gradient;
        gradient = std::move(augmentedGradient);

        MatrixXd augmentedConstraints = MatrixXd::Zero(constraintMatrix.rows(), V);
        augmentedConstraints.leftCols(C) = constraintMatrix;
        constraintMatrix = std::move(augmentedConstraints);

        VectorXd augmentedSeed = VectorXd::Zero(V);
        augmentedSeed.head(C) = seed;
        if(context.previousWarmStart.size() == V)
        {
            augmentedSeed.tail(N) = context.previousWarmStart.tail(N).cwiseMax(0.0);
        }
        seed = std::move(augmentedSeed);
    }

    std::vector<RowVectorXd> forceMaps(static_cast<size_t>(N), RowVectorXd::Zero(V));
    std::vector<double> forceConstants(static_cast<size_t>(N), 0.0);
    for(int stage = 0; stage < N; ++stage)
    {
        MatrixXd positionMap = MatrixXd::Zero(3, C);
        for(int input = 0; input <= stage; ++input)
        {
            positionMap.block<3,3>(
                0, input * context.stageControlDimension) =
                context.dt * context.predictionRotation;
        }

        const Eigen::Matrix4d &boardStage =
            context.parentTransforms[static_cast<size_t>(stage + 1)];
        const Eigen::Vector3d normalStage = detail::inward_contact_normal(
            boardStage.block<3,3>(0,0), _contactParameters.normalAxisInBoard);
        const detail::AffineNormalForce force = detail::affine_normal_force(
            _contactState.measuredNormalForce,
            _contactParameters.forceResponseGain,
            _currentSignedNormalCoordinate,
            normalStage,
            _currentToolPosition,
            boardStage.block<3,1>(0,3),
            positionMap);
        forceMaps[static_cast<size_t>(stage)].head(C) = force.map;
        forceConstants[static_cast<size_t>(stage)] = force.constant;

        if(stage == 0)
        {
            _firstForceMap = force.map;
            _firstForceConstant = force.constant;
            _firstStageNormal = normalStage;
            _firstStageBoardPosition = boardStage.block<3,1>(0,3);
        }

        if(_contactParameters.mode == ContactMode::Loss)
        {
            const double residual = force.constant - _contactParameters.targetForce;
            hessian += 2.0 * _contactParameters.forceWeight
                * forceMaps[static_cast<size_t>(stage)].transpose()
                * forceMaps[static_cast<size_t>(stage)];
            gradient += 2.0 * _contactParameters.forceWeight
                * forceMaps[static_cast<size_t>(stage)].transpose() * residual;
        }
        else
        {
            hessian(C + stage, C + stage) += 2.0 * _contactParameters.slackWeight;
        }
    }

    if(!constraintMode)
    {
        return;
    }

    const int oldRows = constraintMatrix.rows();
    MatrixXd augmentedConstraints = MatrixXd::Zero(oldRows + 4 * N, V);
    augmentedConstraints.topRows(oldRows) = constraintMatrix;
    VectorXd augmentedBounds = VectorXd::Zero(oldRows + 4 * N);
    augmentedBounds.head(oldRows) = constraintVector;

    for(int stage = 0; stage < N; ++stage)
    {
        const int row = oldRows + 4 * stage;
        const int slack = C + stage;
        const RowVectorXd &forceMap = forceMaps[static_cast<size_t>(stage)];
        const double forceConstant = forceConstants[static_cast<size_t>(stage)];

        augmentedConstraints.row(row) = forceMap;
        augmentedConstraints(row, slack) = -1.0;
        augmentedBounds(row) = _contactParameters.targetForce
            + _contactParameters.forceTolerance - forceConstant;

        augmentedConstraints.row(row + 1) = -forceMap;
        augmentedConstraints(row + 1, slack) = -1.0;
        augmentedBounds(row + 1) = forceConstant
            - (_contactParameters.targetForce - _contactParameters.forceTolerance);

        augmentedConstraints(row + 2, slack) = -1.0;
        augmentedBounds(row + 2) = 0.0;
        augmentedConstraints(row + 3, slack) = 1.0;
        augmentedBounds(row + 3) = _contactParameters.maxForceSlack;

        const double seededForce = forceConstant + forceMap.head(C).dot(seed.head(C));
        const double requiredSlack = std::max({
            0.0,
            seededForce - (_contactParameters.targetForce + _contactParameters.forceTolerance),
            (_contactParameters.targetForce - _contactParameters.forceTolerance) - seededForce});
        seed(slack) = std::clamp(
            std::max(seed(slack), requiredSlack), 0.0, _contactParameters.maxForceSlack);
    }

    constraintMatrix = std::move(augmentedConstraints);
    constraintVector = std::move(augmentedBounds);
}

void
SerialLinkContactMPCC::shift_extension_warm_start(
    const MpccQpExtensionContext &context,
    const Eigen::VectorXd &optimum,
    Eigen::VectorXd &shiftedWarmStart)
{
    if(_contactParameters.mode != ContactMode::Constraint || !contact_active())
    {
        return;
    }
    const int C = context.baseControlDimension;
    const int N = context.horizon;
    if(optimum.size() != C + N || shiftedWarmStart.size() != C + N)
    {
        return;
    }
    for(int stage = 0; stage < N - 1; ++stage)
    {
        shiftedWarmStart(C + stage) = optimum(C + stage + 1);
    }
}

void
SerialLinkContactMPCC::on_extended_qp_solution(
    const MpccQpExtensionContext &context,
    const Eigen::VectorXd &optimum)
{
    if(!contact_active())
    {
        return;
    }
    _contactDiagnostics.solverSucceeded = true;
    _contactDiagnostics.fallbackUsed = false;
    _contactDiagnostics.firstCommand.head<3>() =
        context.predictionRotation * optimum.head<3>();
    _contactDiagnostics.firstCommand.tail<3>() =
        context.predictionRotation * optimum.segment<3>(3);
    _contactDiagnostics.predictedFirstStepForce =
        _firstForceConstant + _firstForceMap.dot(optimum.head(context.baseControlDimension));
    if(_contactParameters.mode == ContactMode::Constraint)
    {
        _contactDiagnostics.maxForceSlack =
            std::max(0.0, optimum.tail(context.horizon).maxCoeff());
    }
}

void
SerialLinkContactMPCC::on_twist_resolved(
    const Eigen::Vector<double,6> &commandedBaseTwist,
    const Eigen::Vector<double,6> &realizedBaseTwist,
    const double dt)
{
    if(!contact_active())
    {
        return;
    }
    _contactDiagnostics.firstCommand = commandedBaseTwist;
    _contactDiagnostics.realizedCommand = realizedBaseTwist;
    _contactDiagnostics.twistRealizationError =
        (realizedBaseTwist - commandedBaseTwist).norm();
    const Eigen::Vector3d realizedPosition =
        _currentToolPosition + dt * realizedBaseTwist.head<3>();
    _contactDiagnostics.realizedFirstStepForce =
        _contactState.measuredNormalForce
        + _contactParameters.forceResponseGain
            * (detail::signed_normal_coordinate(
                   _firstStageNormal, realizedPosition, _firstStageBoardPosition)
               - _currentSignedNormalCoordinate);
}

} } // namespace RobotLibrary::Control
