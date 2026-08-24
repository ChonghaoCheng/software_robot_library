/** @file PredictiveContactModel.cpp */

#include <Control/Contact/PredictiveContactModel.h>

#include <Math/MathFunctions.h>

#include <cmath>
#include <stdexcept>

namespace RobotLibrary { namespace Control {
namespace {

Eigen::Matrix3d skew(const Eigen::Vector3d &vector)
{
    Eigen::Matrix3d result;
    result << 0.0, -vector.z(), vector.y(),
              vector.z(), 0.0, -vector.x(),
              -vector.y(), vector.x(), 0.0;
    return result;
}

void validate(const PredictiveContactKinematics &kinematics)
{
    if(kinematics.horizon <= 0 || kinematics.stageControlDimension != 7)
    {
        throw std::invalid_argument(
            "Predictive contact model requires a positive horizon and exactly 7 inputs per stage.");
    }
    if(not std::isfinite(kinematics.dt) || kinematics.dt <= 0.0
       || not std::isfinite(kinematics.expectedDt) || kinematics.expectedDt <= 0.0
       || std::abs(kinematics.dt - kinematics.expectedDt) > 1e-12)
    {
        throw std::invalid_argument(
            "Predictive contact model timestep does not match its identified timestep.");
    }
    if(static_cast<int>(kinematics.parentTransforms.size()) != kinematics.horizon + 1
       || not kinematics.predictionRotation.allFinite()
       || not kinematics.endpointRotationBase.allFinite()
       || not kinematics.endpointPositionBase.allFinite()
       || not kinematics.contactOffsetEndpoint.allFinite()
       || not kinematics.compressionDirectionParent.allFinite()
       || kinematics.compressionDirectionParent.norm() <= 1e-12)
    {
        throw std::invalid_argument("Invalid predictive contact kinematics.");
    }
    for(const Eigen::Matrix4d &transform : kinematics.parentTransforms)
    {
        if(not transform.allFinite())
        {
            throw std::invalid_argument("Predictive parent transforms must be finite.");
        }
    }
    if(kinematics.actuationAware
       && (not std::isfinite(kinematics.realizationAutoregressive)
           || kinematics.realizationAutoregressive < 0.0
           || kinematics.realizationAutoregressive >= 1.0
           || not std::isfinite(kinematics.realizationInputGain)
           || kinematics.realizationInputGain < 0.0
           || kinematics.realizationDelay < 0
           || not std::isfinite(kinematics.initialRealizedRobotNormalVelocity)
           || static_cast<int>(kinematics.pastRobotNormalCommands.size())
                != kinematics.realizationDelay))
    {
        throw std::invalid_argument("Invalid scalar normal realization model/history.");
    }
    for(const double command : kinematics.pastRobotNormalCommands)
    {
        if(not std::isfinite(command))
            throw std::invalid_argument("Past normal commands must be finite.");
    }
}

} // namespace

Eigen::Vector3d
wrist_twist_to_contact_point_velocity(
    const Eigen::Vector<double,6> &wristTwistBase,
    const Eigen::Matrix3d &endpointRotationBase,
    const Eigen::Vector3d &contactOffsetEndpoint)
{
    if(not wristTwistBase.allFinite() || not endpointRotationBase.allFinite()
       || not contactOffsetEndpoint.allFinite())
    {
        throw std::invalid_argument("Wrist/contact-point velocity inputs must be finite.");
    }
    const Eigen::Vector3d offsetBase =
        endpointRotationBase * contactOffsetEndpoint;
    return wristTwistBase.head<3>()
        + wristTwistBase.tail<3>().cross(offsetBase);
}

PredictiveContactAffineModel
build_predictive_contact_affine_model(const PredictiveContactKinematics &kinematics)
{
    validate(kinematics);
    const int N = kinematics.horizon;
    const int D = kinematics.stageControlDimension;
    const int controlDimension = N * D;

    PredictiveContactAffineModel model;
    model.dt = kinematics.dt;
    model.horizon = N;
    model.stageControlDimension = D;
    model.initialContactPointBase = kinematics.endpointPositionBase
        + kinematics.endpointRotationBase * kinematics.contactOffsetEndpoint;
    model.pointVelocityMaps.resize(static_cast<size_t>(N));
    model.pointPositionMaps.resize(static_cast<size_t>(N));
    model.boardLinearVelocitiesBase.resize(static_cast<size_t>(N));
    model.boardAngularVelocitiesBase.resize(static_cast<size_t>(N));
    model.normalsBase.resize(static_cast<size_t>(N));
    model.relativeVelocityMap = Eigen::MatrixXd::Zero(N, controlDimension);
    model.relativeVelocityOffset = Eigen::VectorXd::Zero(N);
    model.commandedRobotNormalVelocityMap =
        Eigen::MatrixXd::Zero(N, controlDimension);
    model.realizedRobotNormalVelocityMap =
        Eigen::MatrixXd::Zero(N, controlDimension);
    model.realizedRobotNormalVelocityOffset = Eigen::VectorXd::Zero(N);

    const Eigen::Vector3d contactOffsetBase =
        kinematics.endpointRotationBase * kinematics.contactOffsetEndpoint;
    Eigen::Matrix<double,3,7> pointVelocity = Eigen::Matrix<double,3,7>::Zero();
    pointVelocity.block<3,3>(0,0) = kinematics.predictionRotation;
    pointVelocity.block<3,3>(0,3) =
        -skew(contactOffsetBase) * kinematics.predictionRotation;

    Eigen::MatrixXd pointPositionMap = Eigen::MatrixXd::Zero(3, controlDimension);
    for(int stage = 0; stage < N; ++stage)
    {
        model.pointVelocityMaps[static_cast<size_t>(stage)] = pointVelocity;
        model.pointPositionMaps[static_cast<size_t>(stage)] = pointPositionMap;

        const Eigen::Matrix4d &board =
            kinematics.parentTransforms[static_cast<size_t>(stage)];
        const Eigen::Matrix4d &boardNext =
            kinematics.parentTransforms[static_cast<size_t>(stage + 1)];
        const Eigen::Matrix3d boardRotation = board.block<3,3>(0,0);
        const Eigen::Vector<double,6> boardBodyTwist =
            RobotLibrary::Math::se3_logarithm(board.inverse() * boardNext)
            / kinematics.dt;
        const Eigen::Vector3d boardLinearBase =
            boardRotation * boardBodyTwist.head<3>();
        const Eigen::Vector3d boardAngularBase =
            boardRotation * boardBodyTwist.tail<3>();
        const Eigen::Matrix3d omega = skew(boardAngularBase);
        const Eigen::Vector3d normal = (
            boardRotation * kinematics.compressionDirectionParent).normalized();
        const Eigen::Vector3d boardPosition = board.block<3,1>(0,3);

        model.boardLinearVelocitiesBase[static_cast<size_t>(stage)] = boardLinearBase;
        model.boardAngularVelocitiesBase[static_cast<size_t>(stage)] = boardAngularBase;
        model.normalsBase[static_cast<size_t>(stage)] = normal;

        Eigen::MatrixXd stageVelocityMap = Eigen::MatrixXd::Zero(3, controlDimension);
        stageVelocityMap.block(0, stage * D, 3, D) = pointVelocity;
        model.commandedRobotNormalVelocityMap.row(stage) =
            normal.transpose() * stageVelocityMap;
        const Eigen::RowVectorXd surfaceVelocityMap =
            normal.transpose() * omega * pointPositionMap;
        const double surfaceVelocityOffset = normal.dot(
            boardLinearBase
            + omega * (model.initialContactPointBase - boardPosition));
        model.relativeVelocityMap.row(stage) =
            model.commandedRobotNormalVelocityMap.row(stage)
            - surfaceVelocityMap;
        model.relativeVelocityOffset(stage) = -surfaceVelocityOffset;

        pointPositionMap.block(0, stage * D, 3, D) +=
            kinematics.dt * pointVelocity;
    }

    if(kinematics.actuationAware)
    {
        Eigen::RowVectorXd realizedMap = Eigen::RowVectorXd::Zero(controlDimension);
        double realizedOffset = kinematics.initialRealizedRobotNormalVelocity;
        for(int stage = 0; stage < N; ++stage)
        {
            model.realizedRobotNormalVelocityMap.row(stage) = realizedMap;
            model.realizedRobotNormalVelocityOffset(stage) = realizedOffset;
            const Eigen::RowVectorXd surfaceMap =
                model.commandedRobotNormalVelocityMap.row(stage)
                - model.relativeVelocityMap.row(stage);
            const double surfaceOffset = -model.relativeVelocityOffset(stage);
            model.relativeVelocityMap.row(stage) = realizedMap - surfaceMap;
            model.relativeVelocityOffset(stage) = realizedOffset - surfaceOffset;

            const int commandStage = stage - kinematics.realizationDelay;
            realizedMap *= kinematics.realizationAutoregressive;
            realizedOffset *= kinematics.realizationAutoregressive;
            if(commandStage >= 0)
            {
                realizedMap += kinematics.realizationInputGain
                    * model.commandedRobotNormalVelocityMap.row(commandStage);
            }
            else
            {
                realizedOffset += kinematics.realizationInputGain
                    * kinematics.pastRobotNormalCommands[static_cast<size_t>(
                        kinematics.realizationDelay + commandStage)];
            }
        }
    }
    else
    {
        model.realizedRobotNormalVelocityMap =
            model.commandedRobotNormalVelocityMap;
    }

    model.penetrationIncrementMap = Eigen::MatrixXd::Zero(N, controlDimension);
    model.penetrationIncrementOffset = Eigen::VectorXd::Zero(N);
    Eigen::RowVectorXd accumulatedMap = Eigen::RowVectorXd::Zero(controlDimension);
    double accumulatedOffset = 0.0;
    for(int lead = 0; lead < N; ++lead)
    {
        accumulatedMap += kinematics.dt * model.relativeVelocityMap.row(lead);
        accumulatedOffset += kinematics.dt * model.relativeVelocityOffset(lead);
        model.penetrationIncrementMap.row(lead) = accumulatedMap;
        model.penetrationIncrementOffset(lead) = accumulatedOffset;
    }
    return model;
}

Eigen::VectorXd
predict_contact_force_affine(const PredictiveContactAffineModel &model,
                             const Eigen::VectorXd &decision,
                             const double measuredForce,
                             const double stiffness)
{
    if(decision.size() != model.horizon * model.stageControlDimension
       || not decision.allFinite() || not std::isfinite(measuredForce)
       || not std::isfinite(stiffness) || stiffness <= 0.0)
    {
        throw std::invalid_argument("Invalid affine contact-force prediction inputs.");
    }
    return Eigen::VectorXd::Constant(model.horizon, measuredForce)
        + stiffness * (model.penetrationIncrementMap * decision
                       + model.penetrationIncrementOffset);
}

PredictiveContactRollout
rollout_predictive_contact_explicit(const PredictiveContactKinematics &kinematics,
                                    const Eigen::VectorXd &decision,
                                    const double measuredForce,
                                    const double stiffness)
{
    validate(kinematics);
    if(decision.size() != kinematics.horizon * kinematics.stageControlDimension
       || not decision.allFinite() || not std::isfinite(measuredForce)
       || not std::isfinite(stiffness) || stiffness <= 0.0)
    {
        throw std::invalid_argument("Invalid explicit contact-force rollout inputs.");
    }
    const PredictiveContactAffineModel maps =
        build_predictive_contact_affine_model(kinematics);
    PredictiveContactRollout result;
    result.relativeNormalVelocity = Eigen::VectorXd::Zero(kinematics.horizon);
    result.commandedRobotNormalVelocity = Eigen::VectorXd::Zero(kinematics.horizon);
    result.realizedRobotNormalVelocity = Eigen::VectorXd::Zero(kinematics.horizon);
    result.penetrationIncrement = Eigen::VectorXd::Zero(kinematics.horizon);
    result.force = Eigen::VectorXd::Zero(kinematics.horizon);
    result.contactPointPositionsBase.resize(static_cast<size_t>(kinematics.horizon + 1));
    Eigen::Vector3d point = maps.initialContactPointBase;
    result.contactPointPositionsBase.front() = point;
    double penetration = 0.0;
    double realizedRobotNormalVelocity =
        kinematics.initialRealizedRobotNormalVelocity;
    for(int stage = 0; stage < kinematics.horizon; ++stage)
    {
        const Eigen::Vector<double,7> input =
            decision.segment<7>(stage * kinematics.stageControlDimension);
        const Eigen::Vector3d pointVelocity =
            maps.pointVelocityMaps[static_cast<size_t>(stage)] * input;
        const double commandedRobotNormalVelocity =
            maps.normalsBase[static_cast<size_t>(stage)].dot(pointVelocity);
        const Eigen::Matrix4d &board =
            kinematics.parentTransforms[static_cast<size_t>(stage)];
        const Eigen::Vector3d boardPosition = board.block<3,1>(0,3);
        const Eigen::Vector3d surfaceVelocity =
            maps.boardLinearVelocitiesBase[static_cast<size_t>(stage)]
            + maps.boardAngularVelocitiesBase[static_cast<size_t>(stage)].cross(
                point - boardPosition);
        const double robotNormalVelocity = kinematics.actuationAware
            ? realizedRobotNormalVelocity : commandedRobotNormalVelocity;
        const double relativeVelocity = robotNormalVelocity
            - maps.normalsBase[static_cast<size_t>(stage)].dot(surfaceVelocity);
        penetration += kinematics.dt * relativeVelocity;
        result.relativeNormalVelocity(stage) = relativeVelocity;
        result.commandedRobotNormalVelocity(stage) = commandedRobotNormalVelocity;
        result.realizedRobotNormalVelocity(stage) = robotNormalVelocity;
        result.penetrationIncrement(stage) = penetration;
        result.force(stage) = measuredForce + stiffness * penetration;
        point += kinematics.dt * pointVelocity;
        result.contactPointPositionsBase[static_cast<size_t>(stage + 1)] = point;
        if(kinematics.actuationAware)
        {
            const int commandStage = stage - kinematics.realizationDelay;
            const double delayedCommand = commandStage >= 0
                ? result.commandedRobotNormalVelocity(commandStage)
                : kinematics.pastRobotNormalCommands[static_cast<size_t>(
                    kinematics.realizationDelay + commandStage)];
            realizedRobotNormalVelocity =
                kinematics.realizationAutoregressive * realizedRobotNormalVelocity
                + kinematics.realizationInputGain * delayedCommand;
        }
    }
    return result;
}

void
add_normalized_predictive_force_cost(const Eigen::MatrixXd &forceMap,
                                     const Eigen::VectorXd &forceConstant,
                                     const double targetForce,
                                     const double weight,
                                     Eigen::MatrixXd &hessian,
                                     Eigen::VectorXd &gradient)
{
    if(forceMap.rows() != forceConstant.size() || forceMap.cols() != hessian.cols()
       || hessian.rows() != hessian.cols() || gradient.size() != hessian.rows()
       || not forceMap.allFinite() || not forceConstant.allFinite()
       || not std::isfinite(targetForce) || targetForce <= 0.0
       || not std::isfinite(weight) || weight < 0.0)
    {
        throw std::invalid_argument("Invalid normalized predictive-force cost inputs.");
    }
    const double scale = weight / (targetForce * targetForce);
    const Eigen::VectorXd residual = forceConstant.array() - targetForce;
    hessian += 2.0 * scale * forceMap.transpose() * forceMap;
    gradient += 2.0 * scale * forceMap.transpose() * residual;
}

} } // namespace RobotLibrary::Control
