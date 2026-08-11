/**
 * @file    SerialLinkMPCC.cpp
 * @brief   MPCC controller implementation for serial link robot arms.
 */

#include <Control/SerialLinkMPCC.h>
#include <Math/CondensedMPC.h>
#include <Math/MathFunctions.h>

#include <Eigen/Geometry>

#include <algorithm>
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
  _dt((dt > 0.0) ? dt : 0.005),
  _qpSolver(parameters.qpsolver)
{
    _controlFrequency = 1.0 / _dt;
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
        "Single-pose tracking is not MPCC. Call set_trajectory() and then step(dt, estimatedProgress).");
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
    _vProgressNominal = std::max(_vProgressMin, 1.0 / duration);
    _vProgressMax = std::max(_vProgressMax, _vProgressNominal);
}

Eigen::VectorXd
SerialLinkMPCC::step(const double dt, const double estimatedProgress)
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
    if(not std::isfinite(estimatedProgress))
    {
        throw std::invalid_argument(
            "[ERROR] [SERIAL LINK MPCC] step(): estimatedProgress must be finite.");
    }

    _dt = dt;
    _pathProgress = std::clamp(estimatedProgress, 0.0, 1.0);
    update();

    const RobotLibrary::Model::Pose referencePose = _trajectory.pose_at_progress(_pathProgress);
    const RobotLibrary::Model::Pose currentPose = endpoint_pose();
    const Eigen::Matrix3d referenceRotation = referencePose.quaternion().toRotationMatrix();

    Eigen::Vector<double,ERROR_DIM> error0;
    error0.head<3>() =
        referenceRotation.transpose() * (currentPose.translation() - referencePose.translation());
    const Eigen::Quaterniond relativeOrientation =
        referencePose.quaternion().conjugate() * currentPose.quaternion();
    error0.tail<3>() = quaternion_to_rotation_vector(relativeOrientation);

    const Eigen::Vector<double,NU> optimalControl = solve_mpcc(error0, referenceRotation);
    _uLast = optimalControl;

    Eigen::Vector<double,6> baseTwist;
    baseTwist.head<3>() = referenceRotation * optimalControl.head<3>();
    baseTwist.tail<3>() = referenceRotation * optimalControl.segment<3>(3);
    return resolve_endpoint_twist(baseTwist);
}

Eigen::Vector<double,SerialLinkMPCC::NU>
SerialLinkMPCC::solve_mpcc(const Eigen::Vector<double,ERROR_DIM> &error0,
                           const Eigen::Matrix3d &referenceRotation)
{
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    const int N = static_cast<int>(_horizon);
    const int errorDim = ERROR_DIM * N;
    const int controlDim = NU * N;
    const double remaining = std::max(0.0, 1.0 - _pathProgress);
    const bool relaxProgressMinimum =
        remaining < static_cast<double>(N) * _dt * _vProgressMin;
    const double progressMinimum = relaxProgressMinimum ? 0.0 : _vProgressMin;

    std::vector<Eigen::Vector<double,ERROR_DIM>> pathTangents(static_cast<size_t>(N));
    double progress = _pathProgress;
    for(int stage = 0; stage < N; ++stage)
    {
        pathTangents[static_cast<size_t>(stage)] =
            path_tangent_at_progress(progress, referenceRotation);

        double rate = _vProgressNominal;
        if(_warmStart.size() == controlDim)
        {
            rate = _warmStart(stage * NU + 6);
        }
        progress = std::clamp(progress + _dt * std::clamp(rate, progressMinimum, _vProgressMax),
                              0.0,
                              1.0);
    }

    // Linearised error prediction:
    // e_j = e_0 + dt * sum_{i<=j}(u_i - tau(s_i) * sdot_i).
    MatrixXd errorResponse = MatrixXd::Zero(errorDim, controlDim);
    VectorXd errorOffset = VectorXd::Zero(errorDim);
    MatrixXd errorWeight = MatrixXd::Zero(errorDim, errorDim);
    for(int stage = 0; stage < N; ++stage)
    {
        const int row = stage * ERROR_DIM;
        errorOffset.segment<ERROR_DIM>(row) = error0;

        const Eigen::Vector<double,ERROR_DIM> &tangent =
            pathTangents[static_cast<size_t>(stage)];
        Eigen::Vector3d tangentDirection = tangent.head<3>();
        if(tangentDirection.norm() > 1e-9)
        {
            tangentDirection.normalize();
        }
        else
        {
            tangentDirection = Eigen::Vector3d::UnitX();
        }

        const Eigen::Matrix3d positionWeight =
            _wContour * Eigen::Matrix3d::Identity()
            + (_wLag - _wContour) * tangentDirection * tangentDirection.transpose();
        errorWeight.block<3,3>(row, row) = positionWeight;
        errorWeight.block<3,3>(row + 3, row + 3) =
            _wOrientation * Eigen::Matrix3d::Identity();

        for(int input = 0; input <= stage; ++input)
        {
            const int column = input * NU;
            errorResponse.block<ERROR_DIM,ERROR_DIM>(row, column) +=
                _dt * Eigen::Matrix<double,ERROR_DIM,ERROR_DIM>::Identity();
            errorResponse.block<ERROR_DIM,1>(row, column + 6) +=
                -_dt * pathTangents[static_cast<size_t>(input)];
        }
    }

    Eigen::Matrix<double,NU,NU> inputWeight =
        Eigen::Matrix<double,NU,NU>::Zero();
    inputWeight.diagonal() <<
        _wInputLinear, _wInputLinear, _wInputLinear,
        _wInputAngular, _wInputAngular, _wInputAngular,
        _wInputProgress;
    const MatrixXd inputWeightHorizon =
        RobotLibrary::Math::block_diagonal(inputWeight, static_cast<unsigned int>(N));

    VectorXd nominalControl = VectorXd::Zero(controlDim);
    MatrixXd pathVelocityResponse = MatrixXd::Zero(errorDim, controlDim);
    const double nominalSeedRate =
        std::min(_vProgressNominal,
                 remaining / std::max(static_cast<double>(N) * _dt, 1e-9));
    for(int stage = 0; stage < N; ++stage)
    {
        const int offset = stage * NU;
        const int row = stage * ERROR_DIM;
        nominalControl.segment<ERROR_DIM>(offset) =
            pathTangents[static_cast<size_t>(stage)] * nominalSeedRate;
        nominalControl(offset + 6) = nominalSeedRate;
        pathVelocityResponse.block<ERROR_DIM,ERROR_DIM>(row, offset).setIdentity();
        pathVelocityResponse.block<ERROR_DIM,1>(row, offset + 6) =
            -pathTangents[static_cast<size_t>(stage)];
    }

    MatrixXd difference = MatrixXd::Zero(controlDim, controlDim);
    for(int stage = 0; stage < N; ++stage)
    {
        difference.block<NU,NU>(stage * NU, stage * NU).setIdentity();
        if(stage > 0)
        {
            difference.block<NU,NU>(stage * NU, (stage - 1) * NU) =
                -Eigen::Matrix<double,NU,NU>::Identity();
        }
    }
    const MatrixXd deltaWeight =
        _wDeltaU * MatrixXd::Identity(controlDim, controlDim);
    VectorXd deltaReference = VectorXd::Zero(controlDim);
    deltaReference.head<NU>() = _uLast;

    MatrixXd H =
        errorResponse.transpose() * errorWeight * errorResponse
        + inputWeightHorizon
        + difference.transpose() * deltaWeight * difference
        + _wVelocityTracking * pathVelocityResponse.transpose() * pathVelocityResponse;
    H += 1e-8 * MatrixXd::Identity(controlDim, controlDim);
    H = 0.5 * (H + H.transpose());

    VectorXd f =
        errorResponse.transpose() * errorWeight * errorOffset
        - difference.transpose() * deltaWeight * deltaReference;
    for(int stage = 0; stage < N; ++stage)
    {
        f(stage * NU + 6) -= _qProgressReward * _dt;
    }

    VectorXd lower = VectorXd::Zero(controlDim);
    VectorXd upper = VectorXd::Zero(controlDim);
    for(int stage = 0; stage < N; ++stage)
    {
        const int offset = stage * NU;
        lower.segment<3>(offset).setConstant(-_vMaxLinear);
        upper.segment<3>(offset).setConstant(_vMaxLinear);
        lower.segment<3>(offset + 3).setConstant(-_vMaxAngular);
        upper.segment<3>(offset + 3).setConstant(_vMaxAngular);
        lower(offset + 6) = progressMinimum;
        upper(offset + 6) = _vProgressMax;
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
        constraintMatrix(box.constraintMatrix.rows(), stage * NU + 6) = _dt;
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
            predictedAdvance += _dt * seed(stage * NU + 6);
        }
        if(predictedAdvance > remaining)
        {
            const double feasibleRate =
                (N > 0 and _dt > 0.0) ? remaining / (static_cast<double>(N) * _dt) : 0.0;
            for(int stage = 0; stage < N; ++stage)
            {
                seed(stage * NU + 6) = feasibleRate;
            }
        }
        return seed;
    };

    VectorXd seed = make_feasible(_warmStart);
    VectorXd optimum;
    try
    {
        optimum = _qpSolver.solve(H, f, constraintMatrix, constraintVector, seed);
        if(optimum.size() != controlDim or not optimum.allFinite())
        {
            throw std::runtime_error("QP returned an invalid solution.");
        }
    }
    catch(const std::exception &)
    {
        optimum = seed;
        if(remaining <= 1e-9)
        {
            for(int stage = 0; stage < N; ++stage)
            {
                optimum.segment<ERROR_DIM>(stage * NU).setZero();
            }
        }
    }
    optimum = make_feasible(optimum);

    _warmStart = optimum;
    for(int stage = 0; stage < N - 1; ++stage)
    {
        _warmStart.segment<NU>(stage * NU) = optimum.segment<NU>((stage + 1) * NU);
    }

    Eigen::Vector<double,NU> firstControl = optimum.head<NU>();
    for(int i = 0; i < 3; ++i)
    {
        firstControl(i) = std::clamp(firstControl(i), -_vMaxLinear, _vMaxLinear);
        firstControl(i + 3) =
            std::clamp(firstControl(i + 3), -_vMaxAngular, _vMaxAngular);
    }
    firstControl(6) = std::clamp(firstControl(6), progressMinimum, _vProgressMax);
    return firstControl;
}

Eigen::Vector<double,6>
SerialLinkMPCC::path_tangent_at_progress(
    const double progress,
    const Eigen::Matrix3d &referenceRotation)
{
    const RobotLibrary::Model::Pose pathPose = _trajectory.pose_at_progress(progress);
    const Eigen::Matrix3d relativeRotation =
        referenceRotation.transpose() * pathPose.quaternion().toRotationMatrix();
    const Eigen::Vector<double,6> bodyTangent =
        _trajectory.tangent_at_progress(progress);

    Eigen::Vector<double,6> localTangent;
    localTangent.head<3>() = relativeRotation * bodyTangent.head<3>();
    localTangent.tail<3>() = relativeRotation * bodyTangent.tail<3>();
    return localTangent;
}

} } // namespace
