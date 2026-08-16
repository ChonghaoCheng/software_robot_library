/**
 * @file    SerialLinkGeometricImpedance.h
 * @author  ChonghaoCheng
 * @date    June 2026
 * @version 1.0.0
 *
 * @brief   Geometric (SE(3)) impedance control with force feedforward / tracking.
 *
 * @details This controller renders a Cartesian spring-damper in a geometrically
 *          consistent way on SE(3) (cf. Seo et al., "Geometric Impedance Control
 *          on SE(3) for Robotic Manipulators", arXiv:2211.07945) and adds an
 *          end-effector wrench feedforward so it can be used as a *force*
 *          controller, not just an impedance.
 *
 *          The control law (all quantities in the base/world frame, [force; torque]
 *          and [linear; angular] ordering) is:
 *
 *              f_impedance = K_p * e_pose + K_d * (v_d - v)
 *              f_force     = F_d + K_f * (F_d - F_meas) + K_fi * integral(F_d - F_meas)
 *              tau         = J^T * (f_impedance + f_force) + N * tau_null + g(q)
 *
 *          With K_f = K_fi = 0 this reduces to geometric impedance plus a pure
 *          force feedforward (the minimal "route A"). When a 6D force/torque
 *          sensor is available, set K_f (and optionally K_fi) to close the loop
 *          on contact force.
 *
 * @copyright (c) 2026
 * @license   OSCL - Free for non-commercial open-source use only.
 */

#ifndef SERIAL_LINK_GEOMETRIC_IMPEDANCE_H
#define SERIAL_LINK_GEOMETRIC_IMPEDANCE_H

#include <Control/Core/SerialLinkBase.h>

#include <memory>

namespace RobotLibrary { namespace Control {

/**
 * @brief Geometric impedance control on SE(3) with wrench feedforward / tracking.
 */
class SerialLinkGeometricImpedance : public SerialLinkBase
{
    public:
        /**
         * @brief Constructor.
         * @param model A pointer to a KinematicTree object.
         * @param endpointName The name of the reference frame in the KinematicTree to be controlled.
         * @param parameters Control gains and algorithm parameters. The Cartesian pose/velocity
         *                   gains are reused as the impedance stiffness K_p and damping K_d.
         */
        SerialLinkGeometricImpedance(std::shared_ptr<RobotLibrary::Model::KinematicTree> model,
                                     const std::string &endpointName,
                                     const RobotLibrary::Control::SerialLinkParameters &parameters = SerialLinkParameters());

        /**
         * @brief Track a Cartesian trajectory while applying the desired/feedforward wrench.
         * @param desiredPose The desired endpoint pose.
         * @param desiredVelocity The desired endpoint twist (6D).
         * @param desiredAcceleration The desired endpoint acceleration (6D, unused in the dynamic-free law).
         * @return The joint torques (nx1).
         */
        Eigen::VectorXd
        track_endpoint_trajectory(const RobotLibrary::Model::Pose &desiredPose,
                                  const Eigen::Vector<double,6>   &desiredVelocity,
                                  const Eigen::Vector<double,6>   &desiredAcceleration) override;

        /**
         * @brief Render impedance damping for a desired endpoint twist, plus the feedforward wrench.
         * @param twist The desired linear & angular velocity (6D).
         * @return The required joint torque.
         */
        Eigen::VectorXd
        resolve_endpoint_twist(const Eigen::Vector<double,6> &twist) override;

        /**
         * @brief Not used in this class. Throws if called.
         */
        Eigen::VectorXd
        resolve_endpoint_motion(const Eigen::Vector<double,6> &endpointMotion) override;

        /**
         * @brief Joint-space PD tracking with optional gravity compensation.
         */
        Eigen::VectorXd
        track_joint_trajectory(const Eigen::VectorXd &desiredPosition,
                               const Eigen::VectorXd &desiredVelocity,
                               const Eigen::VectorXd &desiredAcceleration) override;

        /**
         * @brief Map an end-effector wrench to joint torques, adding the null-space task
         *        and (optionally) gravity compensation, then clamp to effort limits.
         * @param wrench A 6x1 wrench (force 3x1, moment 3x1) the end-effector exerts, base frame.
         * @return The joint torque (nx1).
         */
        Eigen::VectorXd
        map_wrench_to_torque(const Eigen::Vector<double,6> &wrench);

        /**
         * @brief Set the feedforward / target wrench the end-effector should exert on the
         *        environment, expressed in the base frame ([force; torque]).
         */
        void
        set_desired_wrench(const Eigen::Vector<double,6> &wrench) { _desiredWrench = wrench; }

        /**
         * @brief Update the measured wrench EXERTED BY the end-effector ON the environment,
         *        expressed in the base frame. Negate the raw sensor reading if it reports the
         *        wrench applied by the environment on the robot. Used only when force gains > 0.
         */
        void
        update_measured_wrench(const Eigen::Vector<double,6> &wrench) { _measuredWrench = wrench; }

        /**
         * @brief Set the proportional (and optional integral) force-tracking gains. Default zero
         *        gives a pure feedforward; positive gains close the loop on the measured wrench.
         * @param proportional 6x6 gain on the wrench error (F_d - F_meas).
         * @param integral 6x6 gain on the integrated wrench error (default zero).
         */
        void
        set_force_gains(const Eigen::Matrix<double,6,6> &proportional,
                        const Eigen::Matrix<double,6,6> &integral = Eigen::Matrix<double,6,6>::Zero());

        /**
         * @brief Enable or disable gravity compensation (requires a valid dynamic model).
         */
        void
        set_gravity_compensation(bool enable) { _gravityCompensation = enable; }

        /**
         * @brief Reset the force-error integrator.
         */
        void
        reset_force_integral() { _forceIntegral.setZero(); }

        /**
         * @brief Set the desired configuration used for redundancy resolution.
         */
        bool
        set_desired_configuration(const Eigen::VectorXd &configuration);

    protected:

        bool _gravityCompensation = true;                                                           ///< Add g(q) to the commanded torque

        Eigen::VectorXd _desiredConfiguration;                                                      ///< For redundancy resolution

        Eigen::Vector<double,6> _desiredWrench  = Eigen::Vector<double,6>::Zero();                   ///< F_d (end-effector on environment, base frame)

        Eigen::Vector<double,6> _measuredWrench = Eigen::Vector<double,6>::Zero();                   ///< F_meas (end-effector on environment, base frame)

        Eigen::Vector<double,6> _forceIntegral  = Eigen::Vector<double,6>::Zero();                   ///< Integrated wrench error

        Eigen::Matrix<double,6,6> _forceGain         = Eigen::Matrix<double,6,6>::Zero();            ///< K_f
        Eigen::Matrix<double,6,6> _forceIntegralGain = Eigen::Matrix<double,6,6>::Zero();            ///< K_fi

        double _forceIntegralLimit = 50.0;                                                          ///< Anti-windup clamp on each integrator element

        /**
         * @brief Compute the wrench feedforward + tracking term (force law).
         */
        Eigen::Vector<double,6>
        force_wrench();

        /**
         * @brief Compute the per-joint effort limits.
         */
        RobotLibrary::Model::Limits
        compute_control_limits(const unsigned int &jointNumber) override;
};

} } // namespace

#endif
