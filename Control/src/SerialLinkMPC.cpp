/**
 * @file    SerialLinkMPC.cpp
 * @brief   Model predictive Cartesian velocity control for a serial link robot arm.
 *
 * 离散模型：x_{k+1} = x_k + dt * u_k
 * 状态：x = [p; e_R]，p 为末端在基座系下的位置，e_R 为四元数计算的姿态误差；
 * 控制：u = [v; ω]，为末端在基座系下的 twist。
 */

#include <Control/SerialLinkMPC.h>

#include <Eigen/Geometry>

 using Eigen::Matrix;
 using Eigen::VectorXd;
 
 namespace RobotLibrary { namespace Control {
 
 ///////////////////////////////////////////////////////////////
// quaternion_orientation_error: quaternion -> 姿态误差向量
///////////////////////////////////////////////////////////////

Eigen::Vector3d
SerialLinkMPC::quaternion_orientation_error(const Eigen::Quaterniond &qCurrent,
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
 
 ///////////////////////////////////////////////////////////////
 // 构造函数
 ///////////////////////////////////////////////////////////////
 
 SerialLinkMPC::SerialLinkMPC(std::shared_ptr<RobotLibrary::Model::KinematicTree> model,
                              const std::string &endpointName,
                              const RobotLibrary::Control::SerialLinkParameters &parameters,
                              unsigned int horizon,
                              double dt)
     : SerialLinkBase(model, endpointName, parameters),
       _horizon(horizon),
       _dt(dt),
       _innerKinematic(std::make_shared<SerialLinkKinematic>(model, endpointName, parameters)),
       _qpSolver(parameters.qpsolver)
 {
     // 覆盖控制频率：外环 MPC 100 Hz
     _controlFrequency = 100.0;
 
     if (_horizon == 0) _horizon = 1;
     if (_dt <= 0.0)    _dt      = 0.01;
 }
 
 ///////////////////////////////////////////////////////////////
 // 解析控制接口：全部委托给内层 kinematic 控制器
 ///////////////////////////////////////////////////////////////
 
 Eigen::VectorXd
 SerialLinkMPC::resolve_endpoint_motion(const Eigen::Vector<double,6> &endpointMotion)
 {
     return _innerKinematic->resolve_endpoint_motion(endpointMotion);
 }
 
 Eigen::VectorXd
 SerialLinkMPC::resolve_endpoint_twist(const Eigen::Vector<double,6> &twist)
 {
     return _innerKinematic->resolve_endpoint_twist(twist);
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
 
    // 1. 当前状态 x0 = [p; e_R]
    RobotLibrary::Model::Pose currentPose = endpoint_pose();

    Eigen::Vector3d p0 = currentPose.translation();
    Eigen::Quaterniond q0(currentPose.rotation());
    Eigen::Quaterniond qRef(desiredPose.rotation());
    Eigen::Vector3d eR = quaternion_orientation_error(q0, qRef);

    Eigen::Matrix<double,6,1> x0;
    x0.head<3>() = p0;
    x0.tail<3>() = eR;

    // 2. 当前参考状态/控制
    Eigen::Vector3d pRef = desiredPose.translation();

    Eigen::Matrix<double,6,1> xRef;
    xRef.head<3>() = pRef;
    xRef.tail<3>().setZero();
 
     Eigen::Matrix<double,6,1> uRef = desiredVelocity;
 
     // 3. MPC 求解得到最优末端 twist
     Eigen::Matrix<double,6,1> uOpt = solveMPC(x0, xRef, uRef);
 
     // 4. 再进行一次简单限幅
     for (int i = 0; i < 3; ++i)
     {
         uOpt[i] = std::clamp(uOpt[i], -_maxLinearSpeed,  _maxLinearSpeed);
     }
     for (int i = 3; i < 6; ++i)
     {
         uOpt[i] = std::clamp(uOpt[i], -_maxAngularSpeed, _maxAngularSpeed);
     }
 
     // 5. 交给内层解析运动学控制器做 twist -> qdot
     return _innerKinematic->resolve_endpoint_twist(uOpt);
 }
 
 ///////////////////////////////////////////////////////////////
 // 关节空间轨迹：直接复用 kinematic 控制器
 ///////////////////////////////////////////////////////////////
 
 Eigen::VectorXd
 SerialLinkMPC::track_joint_trajectory(const Eigen::VectorXd &desiredPosition,
                                       const Eigen::VectorXd &desiredVelocity,
                                       const Eigen::VectorXd &desiredAcceleration)
 {
     return _innerKinematic->track_joint_trajectory(desiredPosition,
                                                    desiredVelocity,
                                                    desiredAcceleration);
 }
 
 ///////////////////////////////////////////////////////////////
 // 关节速度限制：与 SerialLinkKinematic 相同
 ///////////////////////////////////////////////////////////////
 
 RobotLibrary::Model::Limits
 SerialLinkMPC::compute_control_limits(const unsigned int &jointNumber)
 {
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
             "[ERROR] [SERIAL LINK MPC] compute_control_limits(): "
             "Lower limit for joint " + std::to_string(jointNumber) + " is greater than upper limit ("
             + std::to_string(limits.lower) + " > " + std::to_string(limits.upper) + ").");
     }
 
     return limits;
 }
 
 ///////////////////////////////////////////////////////////////
 // solveMPC: 构造并求解有限时域 QP
 ///////////////////////////////////////////////////////////////
 
 Eigen::Matrix<double,6,1>
 SerialLinkMPC::solveMPC(const Eigen::Matrix<double,6,1> &x0,
                         const Eigen::Matrix<double,6,1> &xRef,
                         const Eigen::Matrix<double,6,1> &uRef)
 {
     using namespace Eigen;
 
     const int nx = 6;
     const int nu = 6;
 
     const unsigned int N  = (_horizon == 0) ? 1 : _horizon;
     const double      dt = (_dt > 0.0) ? _dt : 0.01;
 
     const int stateDim   = nx * static_cast<int>(N);   // N 个状态（x1..xN）堆叠
     const int controlDim = nu * static_cast<int>(N);   // N 个控制（u0..u_{N-1}）堆叠
 
     // x_{k+1} = A x_k + B u_k
     Matrix<double,nx,nx> A = Matrix<double,nx,nx>::Identity();
     Matrix<double,nx,nx> B = dt * Matrix<double,nx,nx>::Identity();
 
     // X = Ax * x0 + Bu * U，X 堆叠的是 x1..xN
     MatrixXd Ax(stateDim, nx);
     MatrixXd Bu(stateDim, controlDim);
     Ax.setZero();
     Bu.setZero();
 
     for (unsigned int k = 0; k < N; ++k)
     {
         Matrix<double,nx,nx> Ak = Matrix<double,nx,nx>::Identity();
         for (unsigned int j = 0; j <= k; ++j)
         {
             Ak = A * Ak;
         }
         Ax.block(static_cast<int>(k) * nx, 0, nx, nx) = Ak;
 
         for (unsigned int j = 0; j <= k; ++j)
         {
             Matrix<double,nx,nx> Phi = Matrix<double,nx,nx>::Identity();
             for (unsigned int m = j + 1; m <= k; ++m)
             {
                 Phi = A * Phi;
             }
             Bu.block(static_cast<int>(k) * nx,
                      static_cast<int>(j) * nu,
                      nx,
                      nu) = Phi * B;
         }
     }
 
     // Q, R 权重
     Matrix<double,nx,nx> Q = Matrix<double,nx,nx>::Zero();
     Q(0,0) = Q(1,1) = Q(2,2) = _wPosition;
     Q(3,3) = Q(4,4) = Q(5,5) = _wOrientation;
 
     Matrix<double,nu,nu> R = Matrix<double,nu,nu>::Zero();
     R(0,0) = R(1,1) = R(2,2) = _wLinearVelocity;
     R(3,3) = R(4,4) = R(5,5) = _wAngularVelocity;
 
     // 块对角 Qb, Rb
     MatrixXd Qb = MatrixXd::Zero(stateDim,   stateDim);
     MatrixXd Rb = MatrixXd::Zero(controlDim, controlDim);
 
     for (unsigned int k = 0; k < N; ++k)
     {
         Qb.block(static_cast<int>(k) * nx,
                  static_cast<int>(k) * nx,
                  nx, nx) = Q;
 
         Rb.block(static_cast<int>(k) * nu,
                  static_cast<int>(k) * nu,
                  nu, nu) = R;
     }
 
     // 参考在整个预测窗内保持不变
     Matrix<double,nx,Dynamic> xref(nx, N);
     Matrix<double,nu,Dynamic> uref(nu, N);
 
     for (unsigned int k = 0; k < N; ++k)
     {
         xref.col(static_cast<int>(k)) = xRef;
         uref.col(static_cast<int>(k)) = uRef;
     }
 
     VectorXd xrefVec = Map<const VectorXd>(xref.data(), stateDim);
     VectorXd urefVec = Map<const VectorXd>(uref.data(), controlDim);
 
     // 代价：0.5 * U^T H U + f^T U （忽略常数项）
     MatrixXd H = Bu.transpose() * Qb * Bu + Rb;
     H += 1e-8 * MatrixXd::Identity(controlDim, controlDim); // 数值阻尼
 
     VectorXd f = Bu.transpose() * Qb * (Ax * x0 - xrefVec) - Rb * urefVec;
 
     // 盒约束：u_min <= u <= u_max
     const int M = controlDim;
     MatrixXd Bineq(2 * M, M);
     Bineq.setZero();
     Bineq.topRows(M).setIdentity();
     Bineq.bottomRows(M) = -MatrixXd::Identity(M, M);
 
     VectorXd umin(M), umax(M);
 
     for (unsigned int k = 0; k < N; ++k)
     {
         const int offset = static_cast<int>(k) * nu;
 
         umin.segment(offset, 3).setConstant(-_maxLinearSpeed);
         umax.segment(offset, 3).setConstant( _maxLinearSpeed);
 
         umin.segment(offset + 3, 3).setConstant(-_maxAngularSpeed);
         umax.segment(offset + 3, 3).setConstant( _maxAngularSpeed);
     }
 
     VectorXd zineq(2 * M);
     zineq.head(M) = umax;
     zineq.tail(M) = -umin;
 
     // warm-start
     if (_warmStart.size() != M)
     {
         _warmStart = urefVec;
     }
 
     VectorXd uSeq;
     try
     {
         // 使用 QPSolver 的不等式约束接口：min 0.5 U'HU + f'U, s.t. B U <= z
         uSeq = _qpSolver.solve(H, f, Bineq, zineq, _warmStart);
     }
     catch (const std::exception &)
     {
         uSeq = urefVec; // 若求解失败则退回参考
     }
 
     if (uSeq.size() != M)
     {
         uSeq = urefVec;
     }
 
     _warmStart = uSeq;
 
     // 取第一个控制量 u0
     Eigen::Matrix<double,6,1> u0 = uSeq.head<6>();
     return u0;
 }
 
 }} // namespace
