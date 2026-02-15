/**
 * @file    SerialLinkMPCC.cpp
 * @brief   MPCC controller implementation for serial link robot arms.
 */

#include <Control/SerialLinkMPCC.h>

#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace RobotLibrary { namespace Control {

Eigen::Vector3d
SerialLinkMPCC::quaternion_orientation_error(const Eigen::Quaterniond &qCurrent,
                                             const Eigen::Quaterniond &qRef)
{
    const double eps = 1e-8;
    Eigen::Quaterniond qErr = qRef * qCurrent.conjugate();
    qErr.normalize();

    if (qErr.w() < 0.0)
    {
        qErr.coeffs() *= -1.0;
    }

    const Eigen::Vector3d v = qErr.vec();
    const double vNorm = v.norm();

    if (!std::isfinite(qErr.w()) || !v.allFinite())
    {
        return Eigen::Vector3d::Zero();
    }

    if (vNorm < eps)
    {
        return Eigen::Vector3d::Zero();
    }

    const double angle = 2.0 * std::atan2(vNorm, qErr.w());
    const Eigen::Vector3d axis = v / vNorm;
    return angle * axis;
}

SerialLinkMPCC::SerialLinkMPCC(std::shared_ptr<RobotLibrary::Model::KinematicTree> model,
                               const std::string &endpointName,
                               const RobotLibrary::Control::SerialLinkParameters &parameters,
                               unsigned int horizon,
                               double dt)
: SerialLinkBase(model, endpointName, parameters),
  _horizon((horizon == 0) ? 1 : horizon),
  _dt((dt > 0.0) ? dt : 0.005),
  _innerKinematic(std::make_shared<SerialLinkKinematic>(model, endpointName, parameters)),
  _qpSolver(parameters.qpsolver)
{
    // MPCC outer-loop rate fixed to 200 Hz as requested.
    _controlFrequency = 200.0;
    _dt = 1.0 / _controlFrequency;
}

Eigen::VectorXd
SerialLinkMPCC::resolve_endpoint_motion(const Eigen::Vector<double,6> &endpointMotion)
{
    return _innerKinematic->resolve_endpoint_motion(endpointMotion);
}

Eigen::VectorXd
SerialLinkMPCC::resolve_endpoint_twist(const Eigen::Vector<double,6> &twist)
{
    return _innerKinematic->resolve_endpoint_twist(twist);
}

Eigen::VectorXd
SerialLinkMPCC::track_joint_trajectory(const Eigen::VectorXd &desiredPosition,
                                       const Eigen::VectorXd &desiredVelocity,
                                       const Eigen::VectorXd &desiredAcceleration)
{
    return _innerKinematic->track_joint_trajectory(desiredPosition,
                                                   desiredVelocity,
                                                   desiredAcceleration);
}

Eigen::VectorXd
SerialLinkMPCC::track_endpoint_trajectory(const RobotLibrary::Model::Pose &desiredPose,
                                          const Eigen::Vector<double,6> &desiredVelocity,
                                          const Eigen::Vector<double,6> &desiredAcceleration)
{
    (void)desiredAcceleration;

    update();

    // Build local frame T0 at the current desired pose.
    const Eigen::Vector3d pRef = desiredPose.translation();
    const Eigen::Matrix3d RRef = desiredPose.quaternion().toRotationMatrix();
    const Eigen::Matrix3d RTW = RRef.transpose();
    const Eigen::Vector3d pTW = -RTW * pRef;

    // Current endpoint pose in T0.
    const RobotLibrary::Model::Pose poseCur = endpoint_pose();
    const Eigen::Vector3d pCur = poseCur.translation();
    const Eigen::Quaterniond qCur = poseCur.quaternion();

    const Eigen::Vector3d pTCur = RTW * pCur + pTW;
    const Eigen::Quaterniond qTW(RTW);
    const Eigen::Quaterniond qTCur = qTW * qCur;
    const Eigen::Vector3d rTCur = quaternion_orientation_error(qTCur, Eigen::Quaterniond::Identity());

    // Desired twist in T0.
    const Eigen::Vector3d vRefB = desiredVelocity.head<3>();
    const Eigen::Vector3d wRefB = desiredVelocity.tail<3>();
    const Eigen::Vector3d vRefT = RTW * vRefB;
    const Eigen::Vector3d wRefT = RTW * wRefB;

    // Progress speed and scaling (kept intentionally).
    const double vProgressRef = std::clamp(_vProgressNominal, 0.0, _vProgressMax);
    double alpha = 1.0;
    if (_vProgressNominal > 1e-6)
    {
        alpha = vProgressRef / _vProgressNominal;
    }
    alpha = std::clamp(alpha, 0.0, 2.0);

    Eigen::Vector<double,NU> uRefScaled = Eigen::Vector<double,NU>::Zero();
    uRefScaled.segment<3>(0) = alpha * vRefT;
    uRefScaled.segment<3>(3) = alpha * wRefT;
    uRefScaled(6) = vProgressRef;

    // Current MPCC state: [x y z rx ry rz s].
    Eigen::Vector<double,NX> x0 = Eigen::Vector<double,NX>::Zero();
    x0.segment<3>(0) = pTCur;
    x0.segment<3>(3) = rTCur;
    x0(6) = _sCurrent;

    // Solve MPCC.
    const Eigen::Vector<double,NU> uOpt = solve_mpcc(x0, uRefScaled);

    // Integrate progress.
    _uLast = uOpt;
    _sCurrent = std::max(0.0, _sCurrent + uOpt(6) * _dt);

    // Map optimal local twist back to base frame.
    const Eigen::Vector3d vB = RTW.transpose() * uOpt.segment<3>(0);
    const Eigen::Vector3d wB = RTW.transpose() * uOpt.segment<3>(3);

    Eigen::Vector<double,6> twistB;
    twistB.head<3>() = vB;
    twistB.tail<3>() = wB;

    return _innerKinematic->resolve_endpoint_twist(twistB);
}

Eigen::Vector<double,SerialLinkMPCC::NU>
SerialLinkMPCC::solve_mpcc(const Eigen::Vector<double,NX> &x0,
                           const Eigen::Vector<double,NU> &uRefScaled)
{
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    const int N = static_cast<int>(_horizon);
    const int stateDim = NX * N;
    const int controlDim = NU * N;

    // Linear model x_{k+1} = x_k + dt*u_k.
    Eigen::Matrix<double,NX,NX> A = Eigen::Matrix<double,NX,NX>::Identity();
    Eigen::Matrix<double,NX,NU> B = _dt * Eigen::Matrix<double,NX,NU>::Identity();

    // Condensed prediction matrices: X = Ax*x0 + Bu*U.
    MatrixXd Ax = MatrixXd::Zero(stateDim, NX);
    MatrixXd Bu = MatrixXd::Zero(stateDim, controlDim);

    for (int k = 0; k < N; ++k)
    {
        Ax.block(k * NX, 0, NX, NX) = A;
        for (int j = 0; j <= k; ++j)
        {
            Bu.block(k * NX, j * NU, NX, NU) = B;
        }
    }

    // Build path-aligned position weight using reference velocity direction.
    Eigen::Vector3d tHat = uRefScaled.head<3>();
    if (tHat.norm() > 1e-6) tHat.normalize();
    else                    tHat = Eigen::Vector3d::UnitX();

    const Eigen::Vector3d aux = (std::abs(tHat.z()) < 0.9)
                              ? Eigen::Vector3d::UnitZ()
                              : Eigen::Vector3d::UnitY();
    const Eigen::Vector3d n1 = (aux.cross(tHat)).normalized();
    const Eigen::Vector3d n2 = (tHat.cross(n1)).normalized();

    Eigen::Matrix3d RPath;
    RPath.col(0) = n1;
    RPath.col(1) = n2;
    RPath.col(2) = tHat;

    const Eigen::Matrix3d QPosPath = Eigen::Vector3d(_wContour, _wContour, _wLag).asDiagonal();
    const Eigen::Matrix3d QPos = RPath * QPosPath * RPath.transpose();

    Eigen::Matrix<double,NX,NX> Qk = Eigen::Matrix<double,NX,NX>::Zero();
    Qk.topLeftCorner<3,3>() = QPos;
    Qk(3,3) = _wOrientation;
    Qk(4,4) = _wOrientation;
    Qk(5,5) = _wOrientation;
    // Qk(6,6) left as zero to let progress be shaped by reward/constraints.

    Eigen::Matrix<double,NU,NU> Rk = Eigen::Matrix<double,NU,NU>::Zero();
    Rk(0,0) = _wInputLinear;
    Rk(1,1) = _wInputLinear;
    Rk(2,2) = _wInputLinear;
    Rk(3,3) = _wInputAngular;
    Rk(4,4) = _wInputAngular;
    Rk(5,5) = _wInputAngular;
    Rk(6,6) = _wInputProgress;

    // Block diagonal Q and R.
    MatrixXd Qb = MatrixXd::Zero(stateDim, stateDim);
    MatrixXd Rb = MatrixXd::Zero(controlDim, controlDim);
    for (int i = 0; i < N; ++i)
    {
        Qb.block(i * NX, i * NX, NX, NX) = Qk;
        Rb.block(i * NU, i * NU, NU, NU) = Rk;
    }

    // Reference stack:
    // position/orientation target in local frame is zero,
    // progress target advances with v_progress_ref.
    VectorXd xRefStack = VectorXd::Zero(stateDim);
    for (int k = 0; k < N; ++k)
    {
        xRefStack(k * NX + 6) = _sCurrent + static_cast<double>(k + 1) * uRefScaled(6) * _dt;
    }

    VectorXd uRefStack = VectorXd::Zero(controlDim);
    for (int k = 0; k < N; ++k)
    {
        uRefStack.segment(k * NU, NU) = uRefScaled;
    }

    // Delta-u regularisation.
    MatrixXd Ed = MatrixXd::Zero(controlDim, controlDim);
    for (int i = 0; i < N; ++i)
    {
        Ed.block(i * NU, i * NU, NU, NU).setIdentity();
        if (i > 0)
        {
            Ed.block(i * NU, (i - 1) * NU, NU, NU) = -Eigen::Matrix<double,NU,NU>::Identity();
        }
    }
    MatrixXd Rdu = _wDeltaU * MatrixXd::Identity(controlDim, controlDim);

    // Velocity tracking term only on first 6 dimensions.
    MatrixXd Rv = MatrixXd::Zero(controlDim, controlDim);
    for (int i = 0; i < N; ++i)
    {
        Rv.block(i * NU, i * NU, 6, 6) = _wVelocityTracking * MatrixXd::Identity(6, 6);
    }

    VectorXd X0 = Ax * x0;
    VectorXd e = X0 - xRefStack;

    // Cost: 0.5 U' H U + f' U.
    MatrixXd H = Bu.transpose() * Qb * Bu + Rb + Ed.transpose() * Rdu * Ed + Rv;
    H += 1e-8 * MatrixXd::Identity(controlDim, controlDim);

    VectorXd f = Bu.transpose() * Qb * e - Rv * uRefStack;

    // Delta-u offset around last applied control.
    VectorXd uOffset = VectorXd::Zero(controlDim);
    for (int i = 0; i < N; ++i)
    {
        uOffset.segment(i * NU, NU) = _uLast;
    }
    f -= Ed.transpose() * Rdu * Ed * uOffset;

    // Progress reward: encourage larger v_progress.
    if (_qProgressReward > 0.0)
    {
        for (int i = 0; i < N; ++i)
        {
            f(i * NU + 6) -= _qProgressReward * _dt;
        }
    }

    // Box constraints on U.
    VectorXd uMin = VectorXd::Zero(controlDim);
    VectorXd uMax = VectorXd::Zero(controlDim);
    for (int i = 0; i < N; ++i)
    {
        const int o = i * NU;
        uMin.segment<3>(o) << -_vMaxLinear, -_vMaxLinear, -_vMaxLinear;
        uMax.segment<3>(o) <<  _vMaxLinear,  _vMaxLinear,  _vMaxLinear;
        uMin.segment<3>(o + 3) << -_vMaxAngular, -_vMaxAngular, -_vMaxAngular;
        uMax.segment<3>(o + 3) <<  _vMaxAngular,  _vMaxAngular,  _vMaxAngular;
        uMin(o + 6) = _vProgressMin;
        uMax(o + 6) = _vProgressMax;
    }

    MatrixXd Bineq = MatrixXd::Zero(2 * controlDim, controlDim);
    Bineq.topRows(controlDim).setIdentity();
    Bineq.bottomRows(controlDim) = -MatrixXd::Identity(controlDim, controlDim);

    VectorXd zIneq(2 * controlDim);
    zIneq.head(controlDim) = uMax;
    zIneq.tail(controlDim) = -uMin;

    // Warm-start.
    if (_warmStart.size() != controlDim)
    {
        _warmStart = uRefStack;
    }

    VectorXd Uopt;
    try
    {
        Uopt = _qpSolver.solve(H, f, Bineq, zIneq, _warmStart);
    }
    catch (const std::exception &)
    {
        Uopt = uRefStack;
    }

    if (Uopt.size() != controlDim)
    {
        Uopt = uRefStack;
    }

    _warmStart = Uopt;

    Eigen::Vector<double,NU> u0 = Uopt.head<NU>();

    // Safety clamp on first control.
    for (int i = 0; i < 3; ++i)
    {
        u0(i) = std::clamp(u0(i), -_vMaxLinear, _vMaxLinear);
        u0(i + 3) = std::clamp(u0(i + 3), -_vMaxAngular, _vMaxAngular);
    }
    u0(6) = std::clamp(u0(6), _vProgressMin, _vProgressMax);

    return u0;
}

RobotLibrary::Model::Limits
SerialLinkMPCC::compute_control_limits(const unsigned int &jointNumber)
{
    // Keep limits consistent with SerialLinkKinematic behaviour.
    RobotLibrary::Model::Limits limits;

    double delta = _model->joint_positions()[jointNumber]
                 - _model->link(jointNumber)->joint().position_limits().lower;

    limits.lower = std::max(-delta * _controlFrequency,
                    std::max(-_model->link(jointNumber)->joint().speed_limit(),
                             -2 * std::sqrt(_maxJointAcceleration * delta)));

    delta = _model->link(jointNumber)->joint().position_limits().upper
          - _model->joint_positions()[jointNumber];

    limits.upper = std::min(delta * _controlFrequency,
                    std::min(_model->link(jointNumber)->joint().speed_limit(),
                             2 * std::sqrt(_maxJointAcceleration * delta)));

    if (limits.lower > limits.upper)
    {
        throw std::runtime_error(
            "[ERROR] [SERIAL LINK MPCC] compute_control_limits(): "
            "Lower limit for joint " + std::to_string(jointNumber) + " is greater than upper limit.");
    }

    return limits;
}

} } // namespace
