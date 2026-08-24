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
}

} // namespace

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
        model.relativeVelocityMap.row(stage) = normal.transpose()
            * (stageVelocityMap - omega * pointPositionMap);
        model.relativeVelocityOffset(stage) = -normal.dot(
            boardLinearBase
            + omega * (model.initialContactPointBase - boardPosition));

        pointPositionMap.block(0, stage * D, 3, D) +=
            kinematics.dt * pointVelocity;
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
    result.penetrationIncrement = Eigen::VectorXd::Zero(kinematics.horizon);
    result.force = Eigen::VectorXd::Zero(kinematics.horizon);
    result.contactPointPositionsBase.resize(static_cast<size_t>(kinematics.horizon + 1));
    Eigen::Vector3d point = maps.initialContactPointBase;
    result.contactPointPositionsBase.front() = point;
    double penetration = 0.0;
    for(int stage = 0; stage < kinematics.horizon; ++stage)
    {
        const Eigen::Vector<double,7> input =
            decision.segment<7>(stage * kinematics.stageControlDimension);
        const Eigen::Vector3d pointVelocity =
            maps.pointVelocityMaps[static_cast<size_t>(stage)] * input;
        const Eigen::Matrix4d &board =
            kinematics.parentTransforms[static_cast<size_t>(stage)];
        const Eigen::Vector3d boardPosition = board.block<3,1>(0,3);
        const Eigen::Vector3d surfaceVelocity =
            maps.boardLinearVelocitiesBase[static_cast<size_t>(stage)]
            + maps.boardAngularVelocitiesBase[static_cast<size_t>(stage)].cross(
                point - boardPosition);
        const double relativeVelocity =
            maps.normalsBase[static_cast<size_t>(stage)].dot(
                pointVelocity - surfaceVelocity);
        penetration += kinematics.dt * relativeVelocity;
        result.relativeNormalVelocity(stage) = relativeVelocity;
        result.penetrationIncrement(stage) = penetration;
        result.force(stage) = measuredForce + stiffness * penetration;
        point += kinematics.dt * pointVelocity;
        result.contactPointPositionsBase[static_cast<size_t>(stage + 1)] = point;
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
