/**
 * @file    SerialLinkMovingFrameMPC.h
 * @brief   Time-indexed MPC for Cartesian trajectories attached to a moving frame.
 */

#ifndef SERIAL_LINK_MOVING_FRAME_MPC_H
#define SERIAL_LINK_MOVING_FRAME_MPC_H

#include <Control/SerialLinkMPC.h>

#include <deque>

namespace RobotLibrary { namespace Control {

/**
 * @brief Predicted pose and twist of a moving reference frame.
 *
 * pose is the transform from the moving frame into the robot base frame.
 * twist is the moving-frame origin twist expressed in the robot base frame:
 * [linear velocity; angular velocity].
 */
struct MovingFrameState
{
    RobotLibrary::Model::Pose pose;
    Eigen::Vector<double,6> twist = Eigen::Vector<double,6>::Zero();
};

/**
 * @brief How the contact normal force is folded into the MPC.
 */
enum class ContactMode
{
    Disabled,    ///< No force term; pure board-glued pose tracking (original behaviour).
    Loss,        ///< Soft force tracking: forceWeight * (F_k - targetForce)^2 added to the cost.
    Constraint   ///< Hard force band F_d +/- tol per step, enforced with per-step slack variables.
};

/**
 * @brief Contact-force parameters for SerialLinkMovingFrameMPC.
 *
 * Defaults are ported from the legacy follow_contact_force server's MpcParameters
 * and config/contact_mpc_constraint.yaml. The board surface is modelled as a local
 * linear spring: F = forceResponseGain * (normal penetration). Normal and tangent
 * axes are given in the board frame (the MovingFrameState pose); they default to the
 * board geometry used in simulation (inward normal along -Y, tangent along +X).
 */
struct ContactParameters
{
    ContactMode mode = ContactMode::Disabled;

    double forceResponseGain   = 1500.0;   ///< Local normal stiffness k_c [N/m].
    double targetForce         = 5.0;      ///< Desired normal force [N].
    double forceTolerance      = 1.0;      ///< Half-width of the force band in Constraint mode [N].
    double forceWeight         = 1.0;      ///< Weight on the force-tracking loss (Loss mode).
    double slackWeight         = 1000.0;   ///< Penalty on per-step force slack (Constraint mode).
    double maxForceSlack       = 50.0;     ///< Upper bound on each force slack [N].

    double tangentPositionWeight = 60.0;   ///< Cost weight on in-surface (tangent) position error.
    double normalPositionWeight  = 3.0;    ///< Cost weight on normal position error (kept small).
    double orientationWeight     = 5.0;    ///< Cost weight on orientation error.
    double velocityWeight        = 0.02;   ///< Regularisation on the commanded twist.
    double deltaVelocityWeight   = 0.1;    ///< Regularisation on twist increments (smoothing).

    /// Board-frame axes defining the contact surface. inwardNormal = -(R_board * normalAxis).
    Eigen::Vector3d normalAxisInBoard  = -Eigen::Vector3d::UnitY();
    Eigen::Vector3d tangentAxisInBoard =  Eigen::Vector3d::UnitX();

    bool   useAcceleration   = true;       ///< Use the constant-acceleration board model (else const velocity).
    double maxBoardAcceleration = 5.0;     ///< Clamp on estimated board linear acceleration [m/s^2].
};

/**
 * @brief MPC controller for tracking a CartesianSpline expressed in a moving frame.
 *
 * The inherited set_trajectory() stores a trajectory in the moving frame. Each
 * control step receives a prediction stack for the moving frame, then converts
 * the relative trajectory into a base-frame reference stack before solving MPC.
 */
class SerialLinkMovingFrameMPC : public SerialLinkMPC
{
    public:
        SerialLinkMovingFrameMPC(
            std::shared_ptr<RobotLibrary::Model::KinematicTree> model,
            const std::string &endpointName,
            const RobotLibrary::Control::SerialLinkParameters &parameters = SerialLinkParameters(),
            unsigned int horizon = 20,
            double dt = 0.002);

        /**
         * @brief Track the stored moving-frame trajectory at the given trajectory time.
         *
         * movingFramePrediction should contain N+1 samples:
         *   sample 0: moving frame at the current control time
         *   sample k: moving frame at current time + k * dt
         * If fewer samples are provided, the last sample is reused.
         */
        Eigen::VectorXd
        track_moving_frame_trajectory_at_time(
            const double &time,
            const std::vector<MovingFrameState> &movingFramePrediction);

        /**
         * @brief Track the stored board-frame trajectory using the in-class board predictor.
         *
         * Builds the board-motion prediction internally (constant acceleration) from the
         * samples supplied through update_board_pose(), then solves the (optionally
         * force-aware) MPC. Equivalent to the two-argument overload but the caller does
         * not have to roll out the prediction itself.
         */
        Eigen::VectorXd
        track_moving_frame_trajectory_at_time(const double &time);

        /**
         * @brief Configure the contact-force behaviour. ContactMode::Disabled (default)
         *        leaves the original pose-tracking behaviour untouched.
         */
        void
        set_contact_parameters(const ContactParameters &parameters) { _contactParameters = parameters; }

        /**
         * @brief Supply the latest measured normal force (F/T wrench projected on the board
         *        normal upstream). Required when the contact mode is not Disabled.
         */
        void
        update_measured_normal_force(double force) { _measuredNormalForce = force; }

        /**
         * @brief Push a new board pose sample for the in-class motion predictor.
         * @param boardPose Transform from the board frame into the robot base frame.
         * @param time      Timestamp of the sample [s].
         */
        void
        update_board_pose(const RobotLibrary::Model::Pose &boardPose, double time);

    private:
        static RobotLibrary::Model::Pose
        compose_reference_pose(const MovingFrameState &frame,
                               const RobotLibrary::Model::Pose &relativePose);

        static Eigen::Vector<double,6>
        compose_reference_twist(const MovingFrameState &frame,
                                const RobotLibrary::Model::Pose &relativePose,
                                const Eigen::Vector<double,6> &relativeTwist);

        /**
         * @brief Roll out the board motion over N+1 samples using the current estimate.
         *        Uses a constant-acceleration model when enabled, else constant velocity.
         */
        std::vector<MovingFrameState>
        predict_board_motion(unsigned int steps) const;

        /**
         * @brief Solve the force-aware MPC (surface-frame pose tracking + board-relative
         *        normal-force loss/constraint) and return the first optimal twist.
         * @param x0         Current tool state [position; rotation_vector].
         * @param xRefStack  Board-glued reference states x1..xN.
         * @param uRefStack  Board-glued feedforward twists u0..u_{N-1}.
         * @param boardStack Predicted board states (index 0 = now, index k = k*dt ahead).
         */
        Eigen::Matrix<double,6,1>
        solve_contact_mpc(const Eigen::Matrix<double,6,1> &x0,
                          const std::vector<Eigen::Matrix<double,6,1>> &xRefStack,
                          const std::vector<Eigen::Matrix<double,6,1>> &uRefStack,
                          const std::vector<MovingFrameState> &boardStack);

        ContactParameters _contactParameters;                 ///< Force behaviour (default Disabled).

        double _measuredNormalForce = 0.0;                    ///< Latest measured normal force [N].

        /// Recent board pose samples (front = oldest) used to estimate velocity & acceleration.
        struct BoardSample { RobotLibrary::Model::Pose pose; double time = 0.0; };
        std::deque<BoardSample> _boardSamples;

        Eigen::Vector3d _boardLinearVelocity     = Eigen::Vector3d::Zero();
        Eigen::Vector3d _boardLinearAcceleration = Eigen::Vector3d::Zero();
        Eigen::Vector3d _boardAngularVelocity    = Eigen::Vector3d::Zero();
        bool _boardEstimateValid = false;

        Eigen::VectorXd _contactWarmStart;                    ///< Warm start for the augmented QP.
        Eigen::Vector<double,6> _previousContactCommand = Eigen::Vector<double,6>::Zero();
        bool _hasPreviousContactCommand = false;
};

} } // namespace

#endif
