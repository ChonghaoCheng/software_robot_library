/**
 * @file PredictiveContactModel.h
 * @brief Pure affine M3 contact prediction for a moving rigid surface.
 */

#ifndef PREDICTIVE_CONTACT_MODEL_H
#define PREDICTIVE_CONTACT_MODEL_H

#include <Eigen/Core>

#include <vector>

namespace RobotLibrary { namespace Control {

struct PredictiveContactKinematics
{
    int horizon = 0;
    int stageControlDimension = 7;
    double dt = 0.002;
    double expectedDt = 0.002;
    Eigen::Matrix3d predictionRotation = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d endpointRotationBase = Eigen::Matrix3d::Identity();
    Eigen::Vector3d endpointPositionBase = Eigen::Vector3d::Zero();
    Eigen::Vector3d contactOffsetEndpoint = Eigen::Vector3d::Zero();
    Eigen::Vector3d compressionDirectionParent{0.0, -1.0, 0.0};
    std::vector<Eigen::Matrix4d> parentTransforms;
    /** Optional scalar commanded-to-realized robot-side normal dynamics. */
    bool actuationAware = false;
    double realizationAutoregressive = 0.0;
    double realizationInputGain = 1.0;
    int realizationDelay = 0;
    double initialRealizedRobotNormalVelocity = 0.0;
    /** Chronological r_cmd[-delay], ..., r_cmd[-1]. */
    std::vector<double> pastRobotNormalCommands;
};

struct PredictiveContactAffineModel
{
    double dt = 0.0;
    int horizon = 0;
    int stageControlDimension = 7;
    Eigen::Vector3d initialContactPointBase = Eigen::Vector3d::Zero();
    std::vector<Eigen::Matrix<double,3,7>> pointVelocityMaps;
    std::vector<Eigen::MatrixXd> pointPositionMaps;
    std::vector<Eigen::Vector3d> boardLinearVelocitiesBase;
    std::vector<Eigen::Vector3d> boardAngularVelocitiesBase;
    std::vector<Eigen::Vector3d> normalsBase;
    Eigen::MatrixXd commandedRobotNormalVelocityMap;
    Eigen::MatrixXd realizedRobotNormalVelocityMap;
    Eigen::VectorXd realizedRobotNormalVelocityOffset;
    Eigen::MatrixXd relativeVelocityMap;
    Eigen::VectorXd relativeVelocityOffset;
    Eigen::MatrixXd penetrationIncrementMap;
    Eigen::VectorXd penetrationIncrementOffset;
};

struct PredictiveContactRollout
{
    Eigen::VectorXd relativeNormalVelocity;
    Eigen::VectorXd commandedRobotNormalVelocity;
    Eigen::VectorXd realizedRobotNormalVelocity;
    Eigen::VectorXd penetrationIncrement;
    Eigen::VectorXd force;
    std::vector<Eigen::Vector3d> contactPointPositionsBase;
};

/** Velocity of the fixed physical contact point from a measured wrist twist. */
Eigen::Vector3d
wrist_twist_to_contact_point_velocity(
    const Eigen::Vector<double,6> &wristTwistBase,
    const Eigen::Matrix3d &endpointRotationBase,
    const Eigen::Vector3d &contactOffsetEndpoint);

/** Build v_rel = A_v U + b_v and Delta-rho = A_rho U + b_rho. */
PredictiveContactAffineModel
build_predictive_contact_affine_model(const PredictiveContactKinematics &kinematics);

/** Condensed M3 force prediction initialized from the current measured force. */
Eigen::VectorXd
predict_contact_force_affine(const PredictiveContactAffineModel &model,
                             const Eigen::VectorXd &decision,
                             double measuredForce,
                             double stiffness);

/** Independent stage-by-stage rigid-body evaluation under frozen endpoint orientation. */
PredictiveContactRollout
rollout_predictive_contact_explicit(const PredictiveContactKinematics &kinematics,
                                    const Eigen::VectorXd &decision,
                                    double measuredForce,
                                    double stiffness);

/** Exact normalized force-cost contribution for 0.5 U' H U + g' U. */
void
add_normalized_predictive_force_cost(const Eigen::MatrixXd &forceMap,
                                     const Eigen::VectorXd &forceConstant,
                                     double targetForce,
                                     double weight,
                                     Eigen::MatrixXd &hessian,
                                     Eigen::VectorXd &gradient);

} } // namespace RobotLibrary::Control

#endif
