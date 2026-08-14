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
 *          The MPC runs at the configured sampling time and outputs a desired
 *          endpoint twist. The shared SerialLinkVelocityBase layer maps this
 *          twist to joint velocities using the same cached model state and
 *          Jacobian as the MPC controller.
 *
 *          The optimisation itself uses the QPSolver class from Math/, with a
 *          condensed prediction model x_{k+1} = x_k + dt * u_k over a finite
 *          horizon, subject to box constraints on the Cartesian velocities.
 */

#ifndef SERIAL_LINK_MPC_H
#define SERIAL_LINK_MPC_H

#include <Control/SerialLinkTimeIndexedMPC.h>
#include <Trajectory/CartesianSpline.h>
#include <Trajectory/CartesianTrajectoryFrame.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

namespace RobotLibrary { namespace Control {

/**
 * @brief MPC-based Cartesian velocity controller for serial link arms.
 *
 * The primary entry point is track_endpoint_trajectory(), which computes a
 * desired endpoint twist over a finite prediction horizon and then maps it to
 * joint velocities using the shared resolved-rate velocity layer.
 */
class SerialLinkMPC : public SerialLinkTimeIndexedMPC
{
    public:

        /**
         * @brief Constructor.
         * @param model        Shared pointer to the kinematic/dynamic model.
         * @param endpointName Name of the endpoint frame on the KinematicTree.
         * @param parameters   Control and QP solver parameters.
         * @param horizon      Prediction horizon length (number of steps).
         * @param dt           MPC sampling time [s] (e.g. 0.002 for 500 Hz).
         */
        SerialLinkMPC(std::shared_ptr<RobotLibrary::Model::KinematicTree> model,
                      const std::string &endpointName,
                      const RobotLibrary::Control::SerialLinkParameters &parameters = SerialLinkParameters(),
                      unsigned int horizon = 20,
                      double dt = 0.002);

        /**
         * @brief Compute joint velocities required to track a Cartesian trajectory.
         *
         * This compatibility entry point only receives one sampled reference, so
         * it constructs a local horizon by linearly extrapolating that reference
         * with the supplied desired velocity.
         */
        Eigen::VectorXd
        track_endpoint_trajectory(const RobotLibrary::Model::Pose &desiredPose,
                                  const Eigen::Vector<double,6>   &desiredVelocity,
                                  const Eigen::Vector<double,6>   &desiredAcceleration) override;

        /**
         * @brief Set the full time-indexed Cartesian reference trajectory.
         *
         * Use this with track_endpoint_trajectory_at_time() when the caller wants
         * the MPC horizon to sample the real future reference poses/twists.
         */
        void
        set_trajectory(const RobotLibrary::Trajectory::CartesianSpline &trajectory) override;

        /** Set the current pose/twist of the trajectory parent frame in base. */
        void
        set_trajectory_frame(
            const RobotLibrary::Trajectory::CartesianTrajectoryFrameState &frame) override;

        /**
         * @brief Clear the currently stored trajectory and warm start.
         */
        void
        clear_trajectory() override;

        /**
         * @brief True when a trajectory has been supplied through set_trajectory().
         */
        bool
        has_trajectory() const override { return _trajectorySet; }

        /**
         * @brief Track the stored time-indexed Cartesian trajectory at the given time.
         *
         * The MPC reference stack is sampled at time + k * dt over the prediction
         * horizon, so the optimisation sees the future path instead of a repeated
         * current reference point.
         */
        Eigen::VectorXd
        track_endpoint_trajectory_at_time(const double &time) override;

        /**
         * @brief Match the unconstrained local feedback bandwidth to target gains.
         *
         * Uses exact discrete-integrator LQR stage and Riccati terminal weights.
         * Baseline behaviour is unchanged until this method is called.
         */
        void set_feedback_bandwidth(double positionGain,
                                    double orientationGain);

    protected:

        /// Length of the MPC prediction horizon (number of time steps).
        unsigned int _horizon = 20;

        /// MPC sampling time [s] (should match 1 / _controlFrequency).
        double _dt = 0.002;

        // State and control weights:
        double _wPosition        = 60.0;   ///< Weight on position error (per axis).
        double _wOrientation     = 10.0;    ///< Weight on orientation (rotation vector) error (per axis).
        double _wLinearVelocity  = 0.1;     ///< Weight on linear velocity deviation (per axis).
        double _wAngularVelocity = 0.1;     ///< Weight on angular velocity deviation (per axis).

        Eigen::Matrix<double,6,6> _terminalStateWeight =
            Eigen::Matrix<double,6,6>::Zero();
        bool _useTerminalStateWeight = false;

        // Cartesian velocity limits:
        double _maxLinearSpeed   = 0.5;     ///< Max linear speed [m/s] for each axis.
        double _maxAngularSpeed  = 0.5;     ///< Max angular speed [rad/s] for each axis.

        /// Optional full reference trajectory used by the time-indexed MPC API.
        RobotLibrary::Trajectory::CartesianSpline _trajectory;

        /// Current trajectory-parent state; identity preserves absolute trajectories.
        RobotLibrary::Trajectory::CartesianTrajectoryFrameState _trajectoryFrame;

        /// True after set_trajectory() has been called.
        bool _trajectorySet = false;

        /// Underlying QP solver instance for the MPC problem.
        QPSolver<double> _qpSolver;

        /// Warm-start vector for the QP decision variable (stacked control inputs).
        Eigen::VectorXd _warmStart;

        /**
         * @brief Solve the MPC problem for a single reference point.
         *
         * This compatibility overload fills the horizon with the same state and
         * control references.
         *
         * @param x0   6x1 current state [position; rotation_vector].
         * @param xRef 6x1 reference state.
         * @param uRef 6x1 reference control (twist).
         * @return 6x1 optimal control input for the first step in the horizon.
         */
        Eigen::Matrix<double,6,1>
        solveMPC(const Eigen::Matrix<double,6,1> &x0,
                 const Eigen::Matrix<double,6,1> &xRef,
                 const Eigen::Matrix<double,6,1> &uRef);

        /**
         * @brief Solve the MPC problem for a full horizon reference stack.
         *
         * @param x0        6x1 current state [position; rotation_vector].
         * @param xRefStack Reference states for x1..xN.
         * @param uRefStack Reference controls for u0..uN-1.
         * @return 6x1 optimal control input for the first step in the horizon.
         */
        Eigen::Matrix<double,6,1>
        solveMPC(const Eigen::Matrix<double,6,1> &x0,
                 const std::vector<Eigen::Matrix<double,6,1>> &xRefStack,
                 const std::vector<Eigen::Matrix<double,6,1>> &uRefStack);
};

} } // namespace

#endif
