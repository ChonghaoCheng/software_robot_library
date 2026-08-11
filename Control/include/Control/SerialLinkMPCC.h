/**
 * @file    SerialLinkMPCC.h
 * @brief   MPCC controller for serial link robot arms.
 *
 * @details This class implements a model predictive contouring control (MPCC)
 *          formulation in a local Cartesian frame. The controller stores the
 *          complete reference path, while the caller supplies the externally
 *          estimated closest progress at every control step.
 *
 *          Error:            e = [position_error; rotation_vector_error]
 *          Control:          u = [linear_velocity; angular_velocity]
 *          Virtual control:  progress speed sdot
 *          Prediction:       e_{k+1} ~= e_k + dt * (u_k - tau(s_k) * sdot_k)
 */

#ifndef SERIAL_LINK_MPCC_H
#define SERIAL_LINK_MPCC_H

#include <Control/SerialLinkVelocityBase.h>
#include <Trajectory/CartesianSpline.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <memory>

namespace RobotLibrary { namespace Control {

class SerialLinkMPCC : public SerialLinkVelocityBase
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

        /**
         * @brief Unsupported single-pose interface required by SerialLinkBase.
         * @throws std::logic_error Always. MPCC requires a complete trajectory.
         */
        Eigen::VectorXd
        track_endpoint_trajectory(const RobotLibrary::Model::Pose &desiredPose,
                                  const Eigen::Vector<double,6>   &desiredVelocity,
                                  const Eigen::Vector<double,6>   &desiredAcceleration) override;

        /**
         * @brief Store a reference path and reset progress/warm-start state.
         */
        void
        set_trajectory(const RobotLibrary::Trajectory::CartesianSpline &trajectory);

        /**
         * @brief Run one MPCC step from an externally estimated closest progress.
         * @param dt Control step [s].
         * @param estimatedProgress Current closest path progress in [0,1].
         */
        Eigen::VectorXd
        step(double dt, double estimatedProgress);

        /**
         * @brief Most recent externally supplied path progress.
         */
        double
        path_progress() const { return _pathProgress; }

    protected:
        /**
         * @brief Generate the local path tangent used by the MPCC prediction.
         *
         * The default implementation uses the 6D SE(3) path tangent. Derived
         * controllers may replace only this geometric term while reusing the
         * same MPCC QP, progress constraints, and warm-start logic.
         */
        virtual Eigen::Vector<double,6>
        path_tangent_at_progress(double progress,
                                 const Eigen::Matrix3d &referenceRotation);

        RobotLibrary::Trajectory::CartesianSpline &
        reference_trajectory() { return _trajectory; }

    private:
        static constexpr int ERROR_DIM = 6;
        static constexpr int NU = 7;

        unsigned int _horizon = 20;
        double _dt = 0.002;

        // Position tracking weights in local path frame
        double _wContour = 1000.0;
        double _wLag = 10.0;
        double _wOrientation = 40.0;

        // Control and regularisation weights
        double _wInputLinear = 3e-4;
        double _wInputAngular = 3e-4;
        double _wInputProgress = 5e-4;
        double _wDeltaU = 3e-5;
        double _wVelocityTracking = 2e-4;
        double _qProgressReward = 2e-4;

        // Control bounds
        double _vMaxLinear = 0.2;
        double _vMaxAngular = 0.2;
        double _vProgressNominal = 0.01;
        double _vProgressMax = 0.1;
        double _vProgressMin = 0.001;

        RobotLibrary::Trajectory::CartesianSpline _trajectory;
        bool _trajectorySet = false;

        // The external estimate is authoritative; predicted progress exists only in the QP.
        double _pathProgress = 0.0;
        Eigen::Vector<double,NU> _uLast = Eigen::Vector<double,NU>::Zero();
        Eigen::VectorXd _warmStart;

        QPSolver<double> _qpSolver;

        Eigen::Vector<double,NU>
        solve_mpcc(const Eigen::Vector<double,ERROR_DIM> &error0,
                   const Eigen::Matrix3d &referenceRotation);
};

} } // namespace

#endif
