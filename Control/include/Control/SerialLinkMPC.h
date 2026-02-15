/**
 * @file    SerialLinkMPC.h
 * @author  (adapted by AI)
 * @date    February 2026
 *
 * @brief   Model predictive Cartesian velocity control for serial link robots.
 *
 * @details This controller implements a Cartesian-space velocity MPC in task space.
 *          The state is the absolute endpoint pose expressed as a 6D vector
 *          [position; rotation_vector] in the base frame, and the control input is
 *          the endpoint twist [linear_velocity; angular_velocity] in the same frame.
 *
 *          The MPC runs at a lower rate (e.g. 100 Hz) and outputs a desired
 *          endpoint twist. An internal SerialLinkKinematic controller is then used
 *          to map this twist to joint velocities, so all joint-space limits,
 *          singularity avoidance, and redundancy handling reuse existing logic.
 *
 *          The optimisation itself uses the QPSolver class from Math/, with a
 *          condensed prediction model x_{k+1} = x_k + dt * u_k over a finite
 *          horizon, subject to box constraints on the Cartesian velocities.
 */

#ifndef SERIAL_LINK_MPC_H
#define SERIAL_LINK_MPC_H

#include <Control/SerialLinkBase.h>
#include <Control/SerialLinkKinematic.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace RobotLibrary { namespace Control {

/**
 * @brief MPC-based Cartesian velocity controller for serial link arms.
 *
 * The primary entry point is track_endpoint_trajectory(), which computes a
 * desired endpoint twist over a finite prediction horizon and then maps it to
 * joint velocities using an internal SerialLinkKinematic object.
 */
class SerialLinkMPC : public SerialLinkBase
{
    public:

        /**
         * @brief Constructor.
         * @param model        Shared pointer to the kinematic/dynamic model.
         * @param endpointName Name of the endpoint frame on the KinematicTree.
         * @param parameters   Control and QP solver parameters.
         * @param horizon      Prediction horizon length (number of steps).
         * @param dt           MPC sampling time [s] (e.g. 0.01 for 100 Hz).
         */
        SerialLinkMPC(std::shared_ptr<RobotLibrary::Model::KinematicTree> model,
                      const std::string &endpointName,
                      const RobotLibrary::Control::SerialLinkParameters &parameters = SerialLinkParameters(),
                      unsigned int horizon = 20,
                      double dt = 0.01);

        /**
         * @brief Solve the joint motion required to achieve a desired endpoint motion.
         *
         * For MPC, this method simply delegates to the internal kinematic controller.
         */
        Eigen::VectorXd
        resolve_endpoint_motion(const Eigen::Vector<double,6> &endpointMotion) override;

        /**
         * @brief Compute joint motion required to execute a specified endpoint twist.
         *
         * For MPC, this also delegates to the internal kinematic controller to
         * preserve the resolved-rate behaviour when MPC is not required.
         */
        Eigen::VectorXd
        resolve_endpoint_twist(const Eigen::Vector<double,6> &twist) override;

        /**
         * @brief Compute joint velocities required to track a Cartesian trajectory.
         *
         * This is the main MPC entry point. The function:
         *  - builds the current 6D state [position; rotation_vector],
         *  - builds a 6D reference [position; rotation_vector],
         *  - runs a finite-horizon MPC to compute the optimal endpoint twist, and
         *  - maps that twist to joint velocities using the internal kinematic
         *    controller.
         */
        Eigen::VectorXd
        track_endpoint_trajectory(const RobotLibrary::Model::Pose &desiredPose,
                                  const Eigen::Vector<double,6>   &desiredVelocity,
                                  const Eigen::Vector<double,6>   &desiredAcceleration) override;

        /**
         * @brief Compute joint motion required to track a joint space trajectory.
         *
         * In the initial implementation this delegates directly to the internal
         * kinematic controller.
         */
        Eigen::VectorXd
        track_joint_trajectory(const Eigen::VectorXd &desiredPosition,
                               const Eigen::VectorXd &desiredVelocity,
                               const Eigen::VectorXd &desiredAcceleration) override;

    protected:

        /**
         * @brief Compute instantaneous limits on joint control.
         *
         * The MPC layer itself constrains Cartesian velocities; joint-space limits
         * are handled identically to SerialLinkKinematic.
         */
        RobotLibrary::Model::Limits
        compute_control_limits(const unsigned int &jointNumber) override;

    private:

        /// Length of the MPC prediction horizon (number of time steps).
        unsigned int _horizon = 20;

        /// MPC sampling time [s] (should match 1 / _controlFrequency).
        double _dt = 0.01;

        // State and control weights:
        double _wPosition        = 60.0;   ///< Weight on position error (per axis).
        double _wOrientation     = 10.0;    ///< Weight on orientation (rotation vector) error (per axis).
        double _wLinearVelocity  = 0.1;     ///< Weight on linear velocity deviation (per axis).
        double _wAngularVelocity = 0.1;     ///< Weight on angular velocity deviation (per axis).

        // Cartesian velocity limits:
        double _maxLinearSpeed   = 0.5;     ///< Max linear speed [m/s] for each axis.
        double _maxAngularSpeed  = 0.5;     ///< Max angular speed [rad/s] for each axis.

        /// Internal kinematic controller used to map twist to joint velocities.
        std::shared_ptr<SerialLinkKinematic> _innerKinematic;

        /// Underlying QP solver instance for the MPC problem.
        QPSolver<double> _qpSolver;

        /// Warm-start vector for the QP decision variable (stacked control inputs).
        Eigen::VectorXd _warmStart;

        /**
         * @brief Quaternion-based orientation error from current to reference.
         *
         * @param qCurrent Current endpoint orientation.
         * @param qRef     Desired endpoint orientation.
         * @return 3x1 orientation error vector (angle-axis form) used by MPC.
         */
        static Eigen::Vector3d
        quaternion_orientation_error(const Eigen::Quaterniond &qCurrent,
                                     const Eigen::Quaterniond &qRef);

        /**
         * @brief Solve the MPC problem for a given initial state and reference.
         *
         * @param x0   6x1 current state [position; orientation_error].
         * @param xRef 6x1 reference state over the horizon.
         * @param uRef 6x1 reference control (twist) over the horizon.
         * @return 6x1 optimal control input for the first step in the horizon.
         */
        Eigen::Matrix<double,6,1>
        solveMPC(const Eigen::Matrix<double,6,1> &x0,
                 const Eigen::Matrix<double,6,1> &xRef,
                 const Eigen::Matrix<double,6,1> &uRef);
};

} } // namespace

#endif
