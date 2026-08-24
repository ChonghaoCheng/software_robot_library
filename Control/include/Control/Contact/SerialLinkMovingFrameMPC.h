/**
 * @file    SerialLinkMovingFrameMPC.h
 * @brief   Time-indexed MPC for Cartesian trajectories attached to a moving frame.
 */

#ifndef SERIAL_LINK_MOVING_FRAME_MPC_H
#define SERIAL_LINK_MOVING_FRAME_MPC_H

#include <Control/Contact/ContactDataStructures.h>
#include <Control/TrajectoryTracking/SerialLinkMPC.h>

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
        set_contact_parameters(const ContactParameters &parameters);

        /**
         * @brief Supply the latest measured normal force (F/T wrench projected on the board
         *        inward normal upstream, positive in compression). Required when
         *        the contact mode is not Disabled.
         */
        void
        update_measured_normal_force(double force);

        /** Latest contact-QP status and force prediction. */
        const ContactMpcDiagnostics&
        contact_diagnostics() const { return _contactDiagnostics; }

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
        predict_board_motion(double controlTime, unsigned int steps) const;

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

        ContactMpcDiagnostics _contactDiagnostics;

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
