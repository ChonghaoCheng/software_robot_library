/**
 * @file    SerialLinkRMPCC.h
 * @brief   Riemannian MPCC controller for serial link robot arms.
 *
 * @details This class implements a Riemannian Model Predictive Contouring Control
 *          (RMPCC) formulation directly on SE(3). The controller stores the
 *          complete time-indexed reference and owns the virtual path progress. The prediction
 *          horizon then optimises future progress along the real path curvature.
 *
 *          The spline and current endpoint pose are first expressed in the
 *          active trajectory parent frame F (board or disturbed_board).
 *
 *          State / error:   E = (^F T_ref(s))^{-1} ^F T,
 *                           e = log_SE3(E) in R^6
 *          Control:         u in R^6 current-endpoint body twist + scalar sdot
 *          Decision vector: z = [u_1..u_N, sdot_1..sdot_N] in R^{7N}
 *          Nominal discrete dynamics:
 *              E_{k+1} = T_ref(s_{k+1})^{-1} T_ref(s_k)
 *                          E_k Exp(dt u_k)
 *              s_{k+1} = s_k + dt sdot_k
 *          The warm-start horizon is rolled out by these SE(3) products. One
 *          stage-wise RTI/SQP linearisation then produces the correction QP.
 *
 *          The first 6 dims of the optimal control are rotated by the current
 *          endpoint orientation into base-frame point velocity [p_dot; omega], then turned into
 *          joint velocities by the shared velocity-control layer.
 *
 *          A constant (rigid) workspace disturbance D is supported via
 *          set_disturbance(): T_active(s) = D * T_ref(s). For a rigid D the body
 *          tangent tau(s) is invariant, so D only shifts the initial error e_0.
 */

#ifndef SERIAL_LINK_RMPCC_H
#define SERIAL_LINK_RMPCC_H

#include <Control/Core/SerialLinkVelocityBase.h>
#include <Control/TrajectoryTracking/RmpccTypes.h>
#include <Trajectory/CartesianSpline.h>
#include <Trajectory/CartesianTrajectoryFrame.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <memory>
#include <string>

namespace RobotLibrary { namespace Control {

/**
 * @brief Tuning parameters for the Riemannian MPCC controller.
 *
 * Pure Eigen/POD struct with no ROS dependency. A front-end (e.g. a ROS action
 * server) fills it from configuration and passes it to the constructor, so the
 * controller stays runtime-tunable without hard-coding weights in the library.
 */
struct RmpccParameters
{
    RmpccParameters() = default;

    static constexpr int TWIST_DIM = 6;

    // Horizon / progress
    int    horizonSteps              = 20;       ///< Prediction horizon length N
    double terminalMultiplier        = 1.0;      ///< Extra weight on the last horizon stage
    bool   autoProgressRate          = true;     ///< Derive progressRateRef from trajectory duration
    double progressRateRef           = 0.0;      ///< Nominal progress rate sdot_ref (1/s)
    double autoProgressRateMin       = 0.0;      ///< Clamp for the auto-derived sdot_ref
    double autoProgressRateMax       = 1.15;     ///< Clamp for the auto-derived sdot_ref
    double progressRateMin           = 0.0;      ///< Lower bound on sdot (0 => derive from ref)
    double progressRateMinMultiplier = 1.0;      ///< Derived lower bound = ref * multiplier
    double progressRateMax           = 0.1;      ///< Upper bound on sdot (0 => derive from ref)
    double progressRateMaxMultiplier = 1.0;      ///< progressRateMax = ref * multiplier when not set
    double progressUpperSlack        = 1e-4;     ///< Slack on the first-step terminal overrun constraint
    double progressScheduleSlack     = 0.0;      ///< Slack on the optional wall-clock schedule limit
    double completionTolerance       = 1e-3;     ///< Treat as "at the end" within this remaining progress

    // Numerics
    double tangentStep               = 1e-3;     ///< Finite-difference step for tau(s)
    double rtiFiniteDifferenceStep   = 1e-6;     ///< Perturbation used for stage-wise RTI Jacobians
    double hessianRegularization     = 1e-8;     ///< Tikhonov term added to the QP Hessian

    // Progress shaping
    double progressReward            = 2e-4;     ///< Stage-integral reward encouraging forward progress
    double progressRateWeight        = 0.0;      ///< Quadratic tracking of sdot toward sdot_ref
    double progressRateSmoothWeight  = 3e-5;     ///< Penalises sdot changes between stages

    // Riemannian tracking
    RmpccPredictorGeometry predictorGeometry = RmpccPredictorGeometry::ExactSE3;
    RmpccObjectiveGeometry objectiveGeometry = RmpccObjectiveGeometry::FullScrewSE3;
    RmpccResidualLinearization residualLinearization =
        RmpccResidualLinearization::FullResidualJacobian;
    double lagWeight                 = 10.0;     ///< Front-end convenience only; controller uses lagWeightMatrix
    double rotationCharacteristicLength = 0.0;  ///< Front-end convenience only; controller uses metric exactly as supplied
    double terminalPositionMultiplier = 1.0;    ///< Extra terminal contour weight for translation
    double terminalRotationMultiplier = 1.0;    ///< Extra terminal contour weight for rotation
    double terminalLagMultiplier      = 1.0;    ///< Extra terminal weight for the scalar lag component
    double terminalLagTranslationScale = 1.0;   ///< Scale translation block of terminal lag cost
    double terminalLagRotationScale    = 1.0;   ///< Scale rotation block of terminal lag cost
    double runningContourScale        = 1.0;    ///< Test-only scale on non-terminal contour costs
    double runningLagScale            = 1.0;    ///< Test-only scale on non-terminal lag costs
    double runningLagTranslationScale = 1.0;    ///< Test-only scale on the translation block of running lag cost
    double runningLagRotationScale    = 1.0;    ///< Test-only scale on the rotation block of running lag cost
    double pathVelocityScale          = 1.0;    ///< Test-only scale on all path-velocity residual costs
    RmpccLagGeometry runningLagGeometry = RmpccLagGeometry::FullScrew; ///< Non-terminal contour/lag projection geometry
    Eigen::Matrix<double, TWIST_DIM, TWIST_DIM> metric =
        Eigen::Matrix<double, TWIST_DIM, TWIST_DIM>::Identity();
    Eigen::Matrix<double, TWIST_DIM, TWIST_DIM> contourWeight =
        (Eigen::Vector<double, TWIST_DIM>() << 1000.0, 1000.0, 1000.0, 40.0, 40.0, 40.0).finished().asDiagonal();
    Eigen::Matrix<double, TWIST_DIM, TWIST_DIM> lagWeightMatrix =
        10.0 * Eigen::Matrix<double, TWIST_DIM, TWIST_DIM>::Identity(); ///< Metric on the projected 6D lag vector
    Eigen::Matrix<double, TWIST_DIM, TWIST_DIM> controlWeight =
        3e-4 * Eigen::Matrix<double, TWIST_DIM, TWIST_DIM>::Identity();
    Eigen::Matrix<double, TWIST_DIM, TWIST_DIM> controlRateWeight =
        3e-5 * Eigen::Matrix<double, TWIST_DIM, TWIST_DIM>::Identity();
    Eigen::Matrix<double, TWIST_DIM, TWIST_DIM> pathVelocityWeight =
        2e-4 * Eigen::Matrix<double, TWIST_DIM, TWIST_DIM>::Identity();

    Eigen::Vector3d linearVelocityMax  = Eigen::Vector3d::Constant(0.2);
    Eigen::Vector3d angularVelocityMax = Eigen::Vector3d::Constant(0.2);
};

/**
 * @brief Per-step diagnostics produced by a SerialLinkRMPCC control step.
 */
struct RmpccDiagnostics
{
    Eigen::Vector<double, 6> bodyTwist   = Eigen::Vector<double, 6>::Zero(); ///< Optimal first body twist u_1
    Eigen::Vector<double, 6> se3Error    = Eigen::Vector<double, 6>::Zero(); ///< e_0 = log(T_ref(s)^-1 T)
    double progressRate          = 0.0;  ///< Optimal first progress rate sdot_1
    double referenceProgress     = 0.0;  ///< Progress s_k used to compute this step's error/control
    double pathProgress          = 0.0;  ///< Integrated progress s_{k+1}
    double se3ErrorNorm          = 0.0;  ///< ||e_0|| (mixed units, diagnostic only)
    double translationError      = 0.0;  ///< ||e_0,translation|| at the same internal progress
    double contourError          = 0.0;  ///< Norm of the contour (perpendicular) error component
    double lagError              = 0.0;  ///< Signed lag (along-path) error component
    double rotationError         = 0.0;  ///< Rotation error magnitude (rad)
    double modelPredictionResidual = 0.0;///< ||e_k - e_k predicted by the preceding local model||
    double twistRealizationError = 0.0;  ///< ||u_QP - R_current^T J qdot_command||
    double qpFirstTwistGradientNorm = 0.0;///< ||(H z + f)_u0||; zero for an interior optimum
    double qpStepNorm           = 0.0;  ///< ||z_opt - z_warm||
    double qpSolveTimeSeconds   = 0.0;  ///< Wall time spent in the task-space QP solve
    double effectiveLoopFrequency = 0.0;///< 1/dt for this controller sample
    double predictedNextErrorNorm = 0.0;///< Norm predicted for the next SE(3) error
    double realizedOneStepErrorNorm = 0.0;///< Norm observed at the next sample
    double activeConstraintCount = 0.0;///< Active task-space QP inequalities
    double linearVelocityLimitActive = 0.0; ///< 1 if any first-stage linear bound is active
    double angularVelocityLimitActive = 0.0;///< 1 if any first-stage angular bound is active
    double jointVelocityLimitActive = 0.0;///< 1 if the realized joint command reaches a bound
    double maxAngularComponent = 0.0; ///< max_i |u_omega_i|
    Eigen::Vector<double,6> feedforwardBodyTwist = Eigen::Vector<double,6>::Zero();
    Eigen::Vector<double,6> correctionBodyTwist = Eigen::Vector<double,6>::Zero();
    double pathVelocityLinearResidual = 0.0; ///< ||u_v - (Ad(E^-1) tau sdot)_v||
    double pathVelocityAngularResidual = 0.0;///< ||u_w - (Ad(E^-1) tau sdot)_w||
    double lagErrorNorm           = 0.0;  ///< ||P_l e_0||, unweighted and unsigned
    double runningContourCost     = 0.0;  ///< Predicted non-terminal sum of weighted contour costs
    double runningLagCost         = 0.0;  ///< Predicted non-terminal sum of weighted lag costs
    double runningPathVelocityCost = 0.0; ///< Predicted horizon path-velocity residual cost
    double terminalContourCost    = 0.0;  ///< Predicted terminal contour contribution
    double terminalLagCost        = 0.0;  ///< Predicted terminal lag contribution
    double terminalLagTranslationCost = 0.0; ///< Translational block contribution to terminal lag cost
    double terminalLagRotationCost = 0.0; ///< Rotational block contribution to terminal lag cost
    double lagTranslationErrorNorm = 0.0; ///< ||(P_l e_0)_translation||
    double lagRotationErrorNorm    = 0.0; ///< ||(P_l e_0)_rotation||
    double runningLagTranslationCost = 0.0; ///< Translational block contribution to running lag cost
    double runningLagRotationCost = 0.0; ///< Rotational block contribution to running lag cost
    double referenceLinearSpeed  = 0.0;  ///< |tau_lin| * sdot_ref
    double referenceAngularSpeed = 0.0;  ///< |tau_ang| * sdot_ref
    double qpStatus              = 1.0;  ///< 1 for a completed step; solver failures throw
    bool   fallbackUsed          = false;///< Compatibility field; strict RMPCC never executes a fallback
};

/**
 * @brief Riemannian MPCC controller using a stored path and internally integrated progress.
 */
class SerialLinkRMPCC : public SerialLinkVelocityBase
{
    public:
        /**
         * @brief Constructor.
         * @param model        Kinematic tree model.
         * @param endpointName Name of the endpoint frame.
         * @param parameters   Shared SerialLinkBase parameters (gains, QP solver options).
         * @param rmpcc        RMPCC-specific tuning parameters.
         */
        SerialLinkRMPCC(std::shared_ptr<RobotLibrary::Model::KinematicTree> model,
                        const std::string &endpointName,
                        const RobotLibrary::Control::SerialLinkParameters &parameters = SerialLinkParameters(),
                        const RmpccParameters &rmpcc = RmpccParameters());

        /**
         * @brief Unsupported single-pose interface required by SerialLinkBase.
         * @throws std::logic_error Always. RMPCC requires a complete trajectory
         *         and an internally integrated progress.
         */
        Eigen::VectorXd
        track_endpoint_trajectory(const RobotLibrary::Model::Pose &desiredPose,
                                  const Eigen::Vector<double,6>   &desiredVelocity,
                                  const Eigen::Vector<double,6>   &desiredAcceleration) override;

        /**
         * @brief Set the reference path and reset progress/warm-start state.
         *        When autoProgressRate is enabled, derives sdot_ref from the duration.
         */
        void
        set_trajectory(const RobotLibrary::Trajectory::CartesianSpline &trajectory);

        /** Set the current rigid trajectory-parent frame in the robot base frame. */
        void
        set_trajectory_frame(
            const RobotLibrary::Trajectory::CartesianTrajectoryFrameState &frame);

        /**
         * @brief Compatibility-only rigid pre-transform D. Default identity.
         *        New integrations should use set_trajectory_frame() and leave D as identity.
         */
        void
        set_disturbance(const Eigen::Matrix4d &disturbance);

        /**
         * @brief Optional wall-clock schedule limit: the first step may not advance
         *        progress beyond this value (+ slack). Pass 1.0 to disable.
         */
        void
        set_schedule_limit(double scheduleProgressLimit);

        /** Lock horizon progress rates to the supplied time-indexed schedule. */
        void
        set_fixed_progress_schedule(bool enabled) { _fixedProgressSchedule = enabled; }

        /**
         * @brief Reset progress, warm start and rate memory (call on a new goal).
         */
        void
        reset();

        /**
         * @brief Most recent internally integrated path progress s in [0, 1].
         */
        double
        path_progress() const { return _pathProgress; }

        /**
         * @brief Reference pose at a given progress, with the current disturbance applied.
         *        Non-const because CartesianSpline::query_state() mutates spline cache state.
         */
        RobotLibrary::Model::Pose
        reference_pose(double progress);

        /** Run one RMPCC step using internally integrated virtual progress. */
        Eigen::VectorXd
        step(double dt);

        /** Compatibility entry point with an externally supplied initial progress. */
        Eigen::VectorXd
        step(double dt, double estimatedProgress);

        /**
         * @brief Diagnostics from the most recent step().
         */
        const RmpccDiagnostics&
        diagnostics() const { return _diagnostics; }

        /** Exact experiment contract for logs and pre-run verification. */
        std::string
        objective_description() const;

    private:
        static constexpr int TWIST_DIM = 6;

        RmpccParameters _rmpcc;
        bool _deriveProgressRateMin = true;
        bool _deriveProgressRateMax = true;

        QPSolver<double> _qpSolver;

        RobotLibrary::Trajectory::CartesianSpline _trajectory;
        bool _trajectorySet = false;

        RobotLibrary::Trajectory::CartesianTrajectoryFrameState _trajectoryFrame;
        Eigen::Matrix4d _disturbance = Eigen::Matrix4d::Identity();
        double _scheduleProgressLimit = 1.0;
        bool _fixedProgressSchedule = false;

        // Controller-owned progress; future progress is predicted inside the QP.
        double _pathProgress = 0.0;
        double _lastProgressRate = 0.0;
        Eigen::Vector<double, TWIST_DIM> _lastBodyTwist = Eigen::Vector<double, TWIST_DIM>::Zero();
        Eigen::Vector<double, TWIST_DIM> _predictedNextError = Eigen::Vector<double, TWIST_DIM>::Zero();
        bool _predictionValid = false;
        Eigen::VectorXd _warmStart;

        RmpccDiagnostics _diagnostics;

        /**
         * @brief Reference transform at a progress value, with disturbance applied.
         *        Nominal path geometry comes from the trajectory; disturbance is applied here.
         *        Non-const: queries the (mutable) spline cache.
         */
        Eigen::Matrix4d
        active_frame_transform_in_base() const;

        /** Reference transform expressed directly in the spline parent frame. */
        Eigen::Matrix4d
        reference_transform_in_trajectory_frame(double progress);

        /** Reference transform expressed in the robot base frame (public reporting only). */
        Eigen::Matrix4d
        reference_transform_in_base(double progress);

        /**
         * @brief Clip a warm-start vector to the box bounds and the first-step progress caps.
         */
        Eigen::VectorXd
        clipped_warm_start(const Eigen::VectorXd &seed,
                           const Eigen::VectorXd &lower,
                           const Eigen::VectorXd &upper,
                           double dt,
                           double remaining,
                           double scheduleRemaining) const;

        /**
         * @brief Roll out the warm-start horizon exactly on SE(3), build one
         *        stage-wise RTI/SQP correction QP, and fill _diagnostics.
         * @param currentTransformInTrajectoryFrame Current endpoint pose in the active
         *        trajectory parent frame as a 4x4 matrix.
         * @param dt               Control step [s].
         */
        void
        solve_rmpcc(const Eigen::Matrix4d &currentTransformInTrajectoryFrame, double dt);
};

} } // namespace

#endif
