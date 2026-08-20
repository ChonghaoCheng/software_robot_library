/**
 * @file    SerialLinkMPC.cpp
 * @brief   Model predictive Cartesian velocity control for a serial link robot arm.
 *
 * 离散模型：x_{k+1} = x_k + dt * u_k
 * 状态：x = [p; r]，p 为末端在基座系下的位置，r 为姿态的 rotation-vector；
 * 控制：u = [v; ω]，为末端在基座系下的 twist。
 */

#include <chrono>
#include <Control/TrajectoryTracking/SerialLinkMPC.h>
#include <Math/CondensedMPC.h>
#include <Math/DiscreteIntegratorLQR.h>

#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <stdexcept>

 using Eigen::Matrix;
 using Eigen::VectorXd;
 using RobotLibrary::Math::quaternion_to_rotation_vector;

namespace RobotLibrary { namespace Control {

namespace {
// Return the shortest relative rotation, avoiding the discontinuity of
// converting two absolute quaternions independently to rotation vectors.
Eigen::Vector3d shortest_rotation_vector(const Eigen::Quaterniond &from,
                                         const Eigen::Quaterniond &to)
{
    Eigen::Quaterniond q = from.conjugate() * to;
    q.normalize();
    if (q.w() < 0.0) q.coeffs() *= -1.0;
    const double w = std::clamp(q.w(), -1.0, 1.0);
    const double angle = 2.0 * std::acos(w);
    const double s = std::sqrt(std::max(0.0, 1.0 - w * w));
    if (s < 1e-8 || angle < 1e-8) return Eigen::Vector3d::Zero();
    return q.vec() * (angle / s);
}
}

 ///////////////////////////////////////////////////////////////
 // 构造函数
 ///////////////////////////////////////////////////////////////
 
 SerialLinkMPC::SerialLinkMPC(std::shared_ptr<RobotLibrary::Model::KinematicTree> model,
                              const std::string &endpointName,
                              const RobotLibrary::Control::SerialLinkParameters &parameters,
                              unsigned int horizon,
                              double dt)
     : SerialLinkTimeIndexedMPC(model, endpointName, parameters),
       _horizon(horizon),
       _dt(dt),
       _qpSolver(parameters.qpsolver)
 {
     if (_horizon == 0) _horizon = 1;
     if (_dt <= 0.0)    _dt      = 0.002;

     _controlFrequency = 1.0 / _dt;
 }

 void
 SerialLinkMPC::set_feedback_bandwidth(const double positionGain,
                                       const double orientationGain)
 {
     _wPosition = RobotLibrary::Math::integrator_stage_weight_for_gain(
         positionGain, _wLinearVelocity, _dt);
     _wOrientation = RobotLibrary::Math::integrator_stage_weight_for_gain(
         orientationGain, _wAngularVelocity, _dt);
     _terminalStateWeight.setZero();
     _terminalStateWeight.diagonal().head<3>().setConstant(
         RobotLibrary::Math::integrator_terminal_weight_for_gain(
             positionGain, _wLinearVelocity, _dt));
     _terminalStateWeight.diagonal().tail<3>().setConstant(
         RobotLibrary::Math::integrator_terminal_weight_for_gain(
             orientationGain, _wAngularVelocity, _dt));
     _useTerminalStateWeight = true;
     _warmStart.resize(0);
 }
 
 ///////////////////////////////////////////////////////////////
 // Time-indexed trajectory API
 ///////////////////////////////////////////////////////////////

 void
 SerialLinkMPC::set_trajectory(const RobotLibrary::Trajectory::CartesianSpline &trajectory)
 {
     _trajectory = trajectory;
     _trajectorySet = true;
     _warmStart.resize(0);
 }

 void
 SerialLinkMPC::set_trajectory_frame(
     const RobotLibrary::Trajectory::CartesianTrajectoryFrameState &frame)
 {
     RobotLibrary::Trajectory::validate_trajectory_frame(frame);
     _trajectoryFrame = frame;
 }

 void
 SerialLinkMPC::clear_trajectory()
 {
     _trajectorySet = false;
     _warmStart.resize(0);
 }

 Eigen::VectorXd
 SerialLinkMPC::track_endpoint_trajectory_at_time(const double &time)
 {
     if (!_trajectorySet)
     {
         throw std::runtime_error(
             "[ERROR] [SERIAL LINK MPC] track_endpoint_trajectory_at_time(): "
             "No trajectory set. Call set_trajectory() first.");
     }

     update();

     const unsigned int N = (_horizon == 0) ? 1 : _horizon;
     const double dt = (_dt > 0.0) ? _dt : 0.002;
     const double startTime = _trajectory.start_time();
     const double endTime = _trajectory.end_time();
     const double t0 = std::isfinite(time) ? std::clamp(time, startTime, endTime) : startTime;

     RobotLibrary::Model::Pose currentPose = endpoint_pose();

     // Keep public diagnostics aligned with the active reference used by this
     // controller, rather than leaving SerialLinkBase's cached error stale.
     const RobotLibrary::Trajectory::CartesianState currentReference =
         RobotLibrary::Trajectory::express_state_in_base(
             _trajectoryFrame, _trajectory.query_state(t0));
     (void)pose_error(currentReference.pose);

     Eigen::Matrix<double,6,1> x0;
     x0.head<3>() = currentPose.translation();
     x0.tail<3>() = quaternion_to_rotation_vector(currentPose.quaternion());

     std::vector<Eigen::Matrix<double,6,1>> xRefStack;
     std::vector<Eigen::Matrix<double,6,1>> uRefStack;
     xRefStack.reserve(N);
     uRefStack.reserve(N);

     for (unsigned int k = 0; k < N; ++k)
     {
         const double stateTime = std::clamp(t0 + static_cast<double>(k + 1) * dt, startTime, endTime);
         const double controlTime = std::clamp(t0 + static_cast<double>(k) * dt, startTime, endTime);

         RobotLibrary::Trajectory::CartesianState stateReference =
             RobotLibrary::Trajectory::express_state_in_base(
                 _trajectoryFrame, _trajectory.query_state(stateTime));
         RobotLibrary::Trajectory::CartesianState controlReference =
             RobotLibrary::Trajectory::express_state_in_base(
                 _trajectoryFrame, _trajectory.query_state(controlTime));

         Eigen::Matrix<double,6,1> xRef;
         xRef.head<3>() = stateReference.pose.translation();
         // Keep the reference on the same local SO(3) branch as the current
         // state. Independent absolute rotation-vector conversions jump at
         // pi (the UR5e ready pose is close to that boundary), producing a
         // spurious 2*pi angular error and driving a joint into its limit.
         const Eigen::Vector3d dtheta_body = shortest_rotation_vector(
             currentPose.quaternion(), stateReference.pose.quaternion());
         xRef.tail<3>() = x0.tail<3>() + currentPose.quaternion().toRotationMatrix() * dtheta_body;

         xRefStack.push_back(xRef);
         uRefStack.push_back(controlReference.twist);
     }

     const auto solveStart = std::chrono::steady_clock::now();
     Eigen::Matrix<double,6,1> uOpt = solveMPC(x0, xRefStack, uRefStack);

     for (int i = 0; i < 3; ++i)
     {
         uOpt[i] = std::clamp(uOpt[i], -_maxLinearSpeed,  _maxLinearSpeed);
     }
     for (int i = 3; i < 6; ++i)
     {
         uOpt[i] = std::clamp(uOpt[i], -_maxAngularSpeed, _maxAngularSpeed);
     }

     // Cartesian MPC clamps directly in base axes, so the clamp-frame twist and
     // the commanded twist are the same vector.
     const Eigen::VectorXd jointCommand = resolve_endpoint_twist(uOpt);
     record_time_indexed_diagnostics(
         uOpt, uOpt, jointCommand, t0, dt, _maxLinearSpeed, _maxAngularSpeed,
         std::chrono::duration<double>(std::chrono::steady_clock::now() - solveStart).count());
     return jointCommand;
 }

 ///////////////////////////////////////////////////////////////
 // track_endpoint_trajectory: MPC 主入口
 ///////////////////////////////////////////////////////////////
 
 Eigen::VectorXd
 SerialLinkMPC::track_endpoint_trajectory(const RobotLibrary::Model::Pose &desiredPose,
                                          const Eigen::Vector<double,6>   &desiredVelocity,
                                          const Eigen::Vector<double,6>   &desiredAcceleration)
 {
     (void)desiredAcceleration; // 速度级控制暂不使用加速度参考
 
     // 更新雅可比与末端姿态
     update();
     (void)pose_error(desiredPose);
 
    // 1. 当前状态 x0 = [p; r]
    RobotLibrary::Model::Pose currentPose = endpoint_pose();

    Eigen::Matrix<double,6,1> x0;
    x0.head<3>() = currentPose.translation();
    x0.tail<3>() = quaternion_to_rotation_vector(currentPose.quaternion());

    // 2. 旧接口只有单个采样点：用当前参考速度外推一个本地 horizon。
    const unsigned int N = (_horizon == 0) ? 1 : _horizon;
    const double dt = (_dt > 0.0) ? _dt : 0.002;
    const Eigen::Vector3d pRef0 = desiredPose.translation();
    const Eigen::Vector3d rRef0 = x0.tail<3>() + currentPose.quaternion().toRotationMatrix() *
        shortest_rotation_vector(currentPose.quaternion(), desiredPose.quaternion());

    std::vector<Eigen::Matrix<double,6,1>> xRefStack;
    std::vector<Eigen::Matrix<double,6,1>> uRefStack;
    xRefStack.reserve(N);
    uRefStack.reserve(N);

    for (unsigned int k = 0; k < N; ++k)
    {
        const double stepTime = static_cast<double>(k + 1) * dt;

        Eigen::Matrix<double,6,1> xRef;
        xRef.head<3>() = pRef0 + stepTime * desiredVelocity.head<3>();
        xRef.tail<3>() = rRef0 + stepTime * desiredVelocity.tail<3>();

        xRefStack.push_back(xRef);
        uRefStack.push_back(desiredVelocity);
    }
 
     // 3. MPC 求解得到最优末端 twist
     Eigen::Matrix<double,6,1> uOpt = solveMPC(x0, xRefStack, uRefStack);
 
     // 4. 再进行一次简单限幅
     for (int i = 0; i < 3; ++i)
     {
         uOpt[i] = std::clamp(uOpt[i], -_maxLinearSpeed,  _maxLinearSpeed);
     }
     for (int i = 3; i < 6; ++i)
     {
         uOpt[i] = std::clamp(uOpt[i], -_maxAngularSpeed, _maxAngularSpeed);
     }
 
     // 5. 使用同一个控制器对象的最新 Jacobian 做 twist -> qdot。
     return resolve_endpoint_twist(uOpt);
 }
 
 ///////////////////////////////////////////////////////////////
 // solveMPC: 构造并求解有限时域 QP
 ///////////////////////////////////////////////////////////////
 
 Eigen::Matrix<double,6,1>
 SerialLinkMPC::solveMPC(const Eigen::Matrix<double,6,1> &x0,
                         const Eigen::Matrix<double,6,1> &xRef,
                         const Eigen::Matrix<double,6,1> &uRef)
 {
     const unsigned int N = (_horizon == 0) ? 1 : _horizon;
     std::vector<Eigen::Matrix<double,6,1>> xRefStack(N, xRef);
     std::vector<Eigen::Matrix<double,6,1>> uRefStack(N, uRef);

     return solveMPC(x0, xRefStack, uRefStack);
 }

 Eigen::Matrix<double,6,1>
 SerialLinkMPC::solveMPC(const Eigen::Matrix<double,6,1> &x0,
                         const std::vector<Eigen::Matrix<double,6,1>> &xRefStack,
                         const std::vector<Eigen::Matrix<double,6,1>> &uRefStack)
 {
     using namespace Eigen;
 
     const int nx = 6;
     const int nu = 6;
 
     const unsigned int N  = (_horizon == 0) ? 1 : _horizon;
     const double      dt = (_dt > 0.0) ? _dt : 0.002;
 
     const int stateDim   = nx * static_cast<int>(N);   // N 个状态（x1..xN）堆叠
     const int controlDim = nu * static_cast<int>(N);   // N 个控制（u0..u_{N-1}）堆叠
 
     // x_{k+1} = A x_k + B u_k ; 凝聚预测 X = Ax * x0 + Bu * U（X 堆叠 x1..xN）
     MatrixXd A = MatrixXd::Identity(nx, nx);
     MatrixXd B = dt * MatrixXd::Identity(nx, nu);

     const RobotLibrary::Math::CondensedPrediction prediction =
         RobotLibrary::Math::condense_prediction(A, B, N);
     const MatrixXd &Ax = prediction.stateTransition;
     const MatrixXd &Bu = prediction.inputResponse;

     // Q, R 权重
     Matrix<double,nx,nx> Q = Matrix<double,nx,nx>::Zero();
     Q(0,0) = Q(1,1) = Q(2,2) = _wPosition;
     Q(3,3) = Q(4,4) = Q(5,5) = _wOrientation;
 
     Matrix<double,nu,nu> R = Matrix<double,nu,nu>::Zero();
     R(0,0) = R(1,1) = R(2,2) = _wLinearVelocity;
     R(3,3) = R(4,4) = R(5,5) = _wAngularVelocity;
 
     // 块对角 Qb, Rb
     MatrixXd Qb = RobotLibrary::Math::block_diagonal(Q, N);
     MatrixXd Rb = RobotLibrary::Math::block_diagonal(R, N);
     if(_useTerminalStateWeight)
     {
         Qb.bottomRightCorner<nx,nx>() = _terminalStateWeight;
     }

     // 参考堆叠：xRefStack 对应 x1..xN，uRefStack 对应 u0..u_{N-1}
     Matrix<double,nx,Dynamic> xref(nx, N);
     Matrix<double,nu,Dynamic> uref(nu, N);

     const Eigen::Matrix<double,6,1> xRefFallback =
         xRefStack.empty() ? x0 : xRefStack.back();
     const Eigen::Matrix<double,6,1> uRefFallback =
         uRefStack.empty() ? Eigen::Matrix<double,6,1>::Zero() : uRefStack.back();
 
     for (unsigned int k = 0; k < N; ++k)
     {
         xref.col(static_cast<int>(k)) =
             (k < xRefStack.size()) ? xRefStack[k] : xRefFallback;
         uref.col(static_cast<int>(k)) =
             (k < uRefStack.size()) ? uRefStack[k] : uRefFallback;
     }
 
     VectorXd xrefVec = Map<const VectorXd>(xref.data(), stateDim);
     VectorXd urefVec = Map<const VectorXd>(uref.data(), controlDim);
 
     // 代价：0.5 * U^T H U + f^T U （忽略常数项）
     MatrixXd H = Bu.transpose() * Qb * Bu + Rb;
     H += 1e-8 * MatrixXd::Identity(controlDim, controlDim); // 数值阻尼
 
     VectorXd f = Bu.transpose() * Qb * (Ax * x0 - xrefVec) - Rb * urefVec;
 
     // 盒约束：u_min <= u <= u_max
     const int M = controlDim;
     VectorXd umin(M), umax(M);

     for (unsigned int k = 0; k < N; ++k)
     {
         const int offset = static_cast<int>(k) * nu;

         umin.segment(offset, 3).setConstant(-_maxLinearSpeed);
         umax.segment(offset, 3).setConstant( _maxLinearSpeed);

         umin.segment(offset + 3, 3).setConstant(-_maxAngularSpeed);
         umax.segment(offset + 3, 3).setConstant( _maxAngularSpeed);
     }

     const RobotLibrary::Math::BoxConstraint box =
         RobotLibrary::Math::box_constraint(umin, umax);
     const MatrixXd &Bineq = box.constraintMatrix;
     const VectorXd &zineq = box.constraintVector;

     // warm-start
     if (_warmStart.size() != M)
     {
         _warmStart = urefVec;
     }
 
     // 使用 QPSolver 的不等式约束接口：min 0.5 U'HU + f'U, s.t. B U <= z
     VectorXd uSeq = _qpSolver.solve(H, f, Bineq, zineq, _warmStart);
     if(uSeq.size() != M or not uSeq.allFinite())
     {
         throw std::runtime_error(
             "[ERROR] [SERIAL LINK MPC] calculate_control(): QP returned an invalid solution.");
     }
     const double constraintViolation = (Bineq * uSeq - zineq).maxCoeff();
     if(not std::isfinite(constraintViolation) or constraintViolation > 1e-6)
     {
         throw std::runtime_error(
             "[ERROR] [SERIAL LINK MPC] calculate_control(): QP returned an infeasible solution (maximum violation "
             + std::to_string(constraintViolation) + ").");
     }
 
     _warmStart = uSeq;
 
     // 取第一个控制量 u0
     Eigen::Matrix<double,6,1> u0 = uSeq.head<6>();
     return u0;
 }
 
 }} // namespace
