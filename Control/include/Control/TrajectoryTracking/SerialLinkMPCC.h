/**
 * @file    SerialLinkMPCC.h
 * @brief   MPCC controller for serial link robot arms.
 *
 * @details This class implements a model predictive contouring control (MPCC)
 *          formulation in a local Cartesian frame. The controller stores the
 *          complete time-indexed reference and owns the virtual path progress.
 *
 *          Error:            e = [position_error; rotation_vector_error]
 *          Control:          u = [linear_velocity; angular_velocity]
 *          Virtual control:  progress speed sdot
 *          Prediction:       e_{k+1} ~= e_k + dt * (u_k - tau(s_k) * sdot_k)
 */

#ifndef SERIAL_LINK_MPCC_H
#define SERIAL_LINK_MPCC_H

#include <Control/Core/SerialLinkVelocityBase.h>
#include <Control/TrajectoryTracking/ParentFrameReferenceMotion.h>
#include <Trajectory/CartesianSpline.h>
#include <Trajectory/CartesianTrajectoryFrame.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <memory>

namespace RobotLibrary { namespace Control {

inline Eigen::Vector<double,6>
mpcc_express_body_tangent_in_prediction_frame(
    const Eigen::Vector<double,6> &bodyTangent,
    const Eigen::Matrix3d &stageRotation,
    const Eigen::Matrix3d &predictionRotation)
{
    const Eigen::Matrix3d rotation =
        predictionRotation.transpose() * stageRotation;
    Eigen::Vector<double,6> result;
    result.head<3>() = rotation * bodyTangent.head<3>();
    result.tail<3>() = rotation * bodyTangent.tail<3>();
    return result;
}

inline Eigen::Matrix3d
mpcc_stage_reference_rotation(
    const Eigen::Matrix4d &predictedParentTransform,
    const Eigen::Matrix4d &pathPose)
{
    return predictedParentTransform.block<3,3>(0,0)
        * pathPose.block<3,3>(0,0);
}

enum class MpccAblationProfile
{
    Baseline,
    FreeProgress,
    OptimizedProgress,
    FeedforwardCorrection,
    MatchedFeedbackGain
};

struct MpccDiagnostics
{
    Eigen::Vector<double,6> error = Eigen::Vector<double,6>::Zero();
    Eigen::Vector<double,6> bodyTwist = Eigen::Vector<double,6>::Zero();
    Eigen::Vector<double,6> realizedBodyTwist = Eigen::Vector<double,6>::Zero();
    double progressRate = 0.0;
    double referenceProgress = 0.0;
    double nextProgress = 0.0;
    double positionError = 0.0;
    double orientationError = 0.0;
    double qpStatus = 1.0;
    double qpIterations = 0.0;
    double qpFinalStepSize = 0.0;
    double qpObjective = 0.0;
    double qpPrimalViolation = 0.0;
    bool qpConverged = false;
    bool qpHitMaxIterations = false;
    double qpActiveSetChanges = 0.0;
    double qpMaximumActiveSetSize = 0.0;
    double qpUniqueActiveConstraints = 0.0;
    double referencePreparationTimeSeconds = 0.0;
    double errorPredictionTimeSeconds = 0.0;
    double costWeightTimeSeconds = 0.0;
    double pathVelocityObjectiveTimeSeconds = 0.0;
    double hessianAssemblyTimeSeconds = 0.0;
    double constraintConstructionTimeSeconds = 0.0;
    double qpSolveTimeSeconds = 0.0;
    double postQpTimeSeconds = 0.0;
    double resolvedRateTimeSeconds = 0.0;
    double totalStepTimeSeconds = 0.0;
    double twistRealizationError = 0.0;
    bool parentFrameMotionActive = false;
    Eigen::Vector<double,6> parentFrameBodyTwist = Eigen::Vector<double,6>::Zero();
    Eigen::Matrix4d measuredParentPose = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d predictedParentPoseFirst = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d predictedParentPoseHorizon = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d parentReferenceFactorFirst = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d repairedReferenceDisplacementFirst = Eigen::Matrix4d::Identity();
    double parentMeasurementTimeSeconds = 0.0;
};

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
                       double dt = 0.0);

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

        /** Configure a reproducible objective/progress ablation before set_trajectory(). */
        void
        set_ablation_profile(MpccAblationProfile profile);

        /** Lock horizon progress rates to the supplied time-indexed schedule. */
        void
        set_fixed_progress_schedule(bool enabled) { _fixedProgressSchedule = enabled; }

        /**
         * @brief Set the Cartesian angular-velocity bound used by the MPCC QP.
         *
         * This is an explicit experiment/profile override. The production
         * default remains 0.2 rad/s.
         */
        void
        set_angular_velocity_limit(double angularVelocityMax);

        /**
         * @brief Set the Cartesian linear-velocity bound used by the MPCC QP.
         *
         * This is an explicit experiment/profile override. The production
         * default remains 0.2 m/s.
         */
        void
        set_linear_velocity_limit(double linearVelocityMax);

        /** Cartesian linear-velocity bound used by the MPCC QP. */
        double
        linear_velocity_limit() const { return _vMaxLinear; }

        /** Cartesian angular-velocity bound used by the MPCC QP. */
        double
        angular_velocity_limit() const { return _vMaxAngular; }

        /** Set the rigid trajectory parent pose/twist in the robot base frame. */
        void
        set_trajectory_frame(
            const RobotLibrary::Trajectory::CartesianTrajectoryFrameState &frame);

        /** Set a timestamped parent-frame measurement for causal prediction. */
        void
        set_trajectory_frame(
            const RobotLibrary::Trajectory::CartesianTrajectoryFrameState &frame,
            double timestampSeconds);

        /** Run one MPCC step using internally integrated virtual progress. */
        Eigen::VectorXd
        step(double dt);

        /** Compatibility entry point with an externally supplied initial progress. */
        Eigen::VectorXd
        step(double dt, double estimatedProgress);

        /**
         * @brief Current internally integrated path progress.
         */
        double
        path_progress() const { return _pathProgress; }

        const MpccDiagnostics &
        diagnostics() const { return _diagnostics; }

        /** Enable detailed steady-clock substages for timing experiments. */
        void set_timing_diagnostics(bool enabled) { _timingDiagnosticsEnabled = enabled; }

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

        const RobotLibrary::Trajectory::CartesianTrajectoryFrameState &
        trajectory_frame() const { return _trajectoryFrame; }

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

        // Local small-error bandwidth of the classic controller.  Its
        // quaternion-vector orientation error is approximately theta/2.
        double _matchedPositionGain = 20.0;
        double _matchedOrientationGain = 2.5;

        MpccAblationProfile _ablationProfile = MpccAblationProfile::Baseline;

        RobotLibrary::Trajectory::CartesianSpline _trajectory;
        RobotLibrary::Trajectory::CartesianTrajectoryFrameState _trajectoryFrame;
        CausalParentFrameMotion _parentFrameMotion;
        bool _trajectorySet = false;

        // Controller-owned progress; future progress is predicted inside the QP.
        double _pathProgress = 0.0;
        bool _fixedProgressSchedule = false;
        Eigen::Vector<double,NU> _uLast = Eigen::Vector<double,NU>::Zero();
        Eigen::VectorXd _warmStart;

        QPSolver<double> _qpSolver;
        SolverOptions<double> _qpOptions;
        MpccDiagnostics _diagnostics;
        bool _timingDiagnosticsEnabled = false;
        std::uint64_t _controlStepIndex = 0;

        Eigen::Vector<double,NU>
        solve_mpcc(const Eigen::Vector<double,ERROR_DIM> &error0,
                   const Eigen::Matrix3d &referenceRotation);
};

} } // namespace

#endif
