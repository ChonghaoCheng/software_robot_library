/**
 * @file    SerialLinkMPCC.h
 * @brief   MPCC controller for serial link robot arms.
 *
 * @details This class implements a model predictive contouring control (MPCC)
 *          formulation in Cartesian space with a 7D state:
 *              x = [position(3), rotation_vector(3), progress]
 *          and a 7D control:
 *              u = [linear_velocity(3), angular_velocity(3), progress_speed].
 *
 *          The optimisation runs in a local path-aligned frame built from the
 *          current desired pose/twist, and the first 6 dimensions of the
 *          resulting control are mapped back to the base frame before being
 *          converted to joint velocities with an internal SerialLinkKinematic.
 */

#ifndef SERIAL_LINK_MPCC_H
#define SERIAL_LINK_MPCC_H

#include <Control/SerialLinkBase.h>
#include <Control/SerialLinkKinematic.h>

#include <Eigen/Core>
#include <memory>

namespace RobotLibrary { namespace Control {

class SerialLinkMPCC : public SerialLinkBase
{
    public:
        /**
         * @brief Constructor.
         * @param model Kinematic tree model.
         * @param endpointName Name of endpoint frame.
         * @param parameters Shared control parameters.
         * @param horizon Number of prediction steps.
         * @param dt Sampling time in seconds.
         */
        SerialLinkMPCC(std::shared_ptr<RobotLibrary::Model::KinematicTree> model,
                       const std::string &endpointName,
                       const RobotLibrary::Control::SerialLinkParameters &parameters = SerialLinkParameters(),
                       unsigned int horizon = 20,
                       double dt = 0.005);

        Eigen::VectorXd
        resolve_endpoint_motion(const Eigen::Vector<double,6> &endpointMotion) override;

        Eigen::VectorXd
        resolve_endpoint_twist(const Eigen::Vector<double,6> &twist) override;

        Eigen::VectorXd
        track_endpoint_trajectory(const RobotLibrary::Model::Pose &desiredPose,
                                  const Eigen::Vector<double,6>   &desiredVelocity,
                                  const Eigen::Vector<double,6>   &desiredAcceleration) override;

        Eigen::VectorXd
        track_joint_trajectory(const Eigen::VectorXd &desiredPosition,
                               const Eigen::VectorXd &desiredVelocity,
                               const Eigen::VectorXd &desiredAcceleration) override;

    protected:
        RobotLibrary::Model::Limits
        compute_control_limits(const unsigned int &jointNumber) override;

    private:
        // MPCC dimensions
        static constexpr int NX = 7;
        static constexpr int NU = 7;

        unsigned int _horizon = 20;
        double _dt = 0.005;

        // Position tracking weights in local path frame
        double _wContour = 150.0;
        double _wLag = 10.0;
        double _wOrientation = 1.0;

        // Control and regularisation weights
        double _wInputLinear = 1e-3;
        double _wInputAngular = 1e-3;
        double _wInputProgress = 1e-2;
        double _wDeltaU = 1e-2;
        double _wVelocityTracking = 5e-3;
        double _qProgressReward = 1e-7;

        // Control bounds
        double _vMaxLinear = 1.0;
        double _vMaxAngular = 0.2;
        double _vProgressNominal = 0.1;
        double _vProgressMax = 0.2;
        double _vProgressMin = 1e-3;

        // MPCC internal state
        double _sCurrent = 0.0;
        Eigen::Vector<double,NU> _uLast = Eigen::Vector<double,NU>::Zero();
        Eigen::VectorXd _warmStart;

        std::shared_ptr<SerialLinkKinematic> _innerKinematic;
        QPSolver<double> _qpSolver;

        static Eigen::Vector3d
        rotation_log(const Eigen::Matrix3d &R);

        Eigen::Vector<double,NU>
        solve_mpcc(const Eigen::Vector<double,NX> &x0,
                   const Eigen::Vector<double,NU> &uRefScaled);
};

} } // namespace

#endif

