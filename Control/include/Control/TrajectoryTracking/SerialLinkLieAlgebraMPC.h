/**
 * @file SerialLinkLieAlgebraMPC.h
 * @brief Time-indexed Cartesian MPC in SE(3) logarithmic error coordinates.
 */

#ifndef SERIAL_LINK_LIE_ALGEBRA_MPC_H
#define SERIAL_LINK_LIE_ALGEBRA_MPC_H

#include <Control/TrajectoryTracking/SerialLinkTimeIndexedMPC.h>

#include <Eigen/Core>
#include <vector>

namespace RobotLibrary { namespace Control {

/**
 * @brief Convex error-state MPC for a fully actuated endpoint twist on SE(3).
 *
 * Unlike SerialLinkMPC, the state is the relative logarithmic pose error
 * Log(T_ref^-1 T), not an absolute [position; rotation-vector] coordinate.
 * Unlike RMPCC, the reference remains time-indexed and there is no optimized
 * path-progress state.
 */
class SerialLinkLieAlgebraMPC : public SerialLinkTimeIndexedMPC
{
    public:
        SerialLinkLieAlgebraMPC(
            std::shared_ptr<RobotLibrary::Model::KinematicTree> model,
            const std::string &endpointName,
            const RobotLibrary::Control::SerialLinkParameters &parameters = SerialLinkParameters(),
            unsigned int horizon = 20,
            double dt = 0.002);

        Eigen::VectorXd
        track_endpoint_trajectory(
            const RobotLibrary::Model::Pose &desiredPose,
            const Eigen::Vector<double,6> &desiredVelocity,
            const Eigen::Vector<double,6> &desiredAcceleration) override;

        void
        set_trajectory(const RobotLibrary::Trajectory::CartesianSpline &trajectory) override;

        void
        set_trajectory_frame(
            const RobotLibrary::Trajectory::CartesianTrajectoryFrameState &frame) override;

        void clear_trajectory() override;
        bool has_trajectory() const override { return _trajectorySet; }

        Eigen::VectorXd
        track_endpoint_trajectory_at_time(const double &time) override;

        /** Configure exact discrete-integrator LQR stage/terminal weights. */
        void set_feedback_bandwidth(double positionGain,
                                    double orientationGain);

    private:
        Eigen::Matrix<double,6,1>
        solve_error_mpc(
            const Eigen::Matrix<double,6,1> &initialError,
            const std::vector<Eigen::Matrix<double,6,1>> &referenceBodyTwists);

        unsigned int _horizon{20};
        double _dt{0.002};

        double _wPosition{60.0};
        double _wOrientation{10.0};
        double _wLinearVelocity{0.1};
        double _wAngularVelocity{0.1};
        Eigen::Matrix<double,6,6> _terminalStateWeight =
            Eigen::Matrix<double,6,6>::Zero();
        bool _useTerminalStateWeight{false};
        double _maxLinearSpeed{0.5};
        double _maxAngularSpeed{0.5};

        RobotLibrary::Trajectory::CartesianSpline _trajectory;
        RobotLibrary::Trajectory::CartesianTrajectoryFrameState _trajectoryFrame;
        bool _trajectorySet{false};

        QPSolver<double> _qpSolver;
        Eigen::VectorXd _warmStart;
};

} } // namespace RobotLibrary::Control

#endif
