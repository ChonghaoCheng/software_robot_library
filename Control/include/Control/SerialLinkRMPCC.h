/**
 * @file    SerialLinkRMPCC.h
 * @brief   Riemannian MPCC controller for serial link robot arms.
 *
 * @details This class implements a Riemannian Model Predictive Contouring Control
 *          (RMPCC) formulation directly on SE(3). The controller stores the
 *          complete reference path, while the caller supplies the externally
 *          estimated closest progress at every control step. The prediction
 *          horizon then optimises future progress along the real path curvature.
 *
 *          State / error:   e = log_SE3( T_ref(s)^{-1} T ) in R^6  (body twist)
 *          Control:         u in R^6 body twist + scalar progress rate sdot
 *          Decision vector: z = [u_1..u_N, sdot_1..sdot_N] in R^{7N}
 *          Prediction:      e_j ~= e_0 + dt * sum_{i<=j} ( u_i - tau(s_i) * sdot_i )
 *
 *          The first 6 dims of the optimal control are rotated into the base
 *          frame as endpoint point velocity [p_dot; omega], then turned into
 *          joint velocities by the shared velocity-control layer.
 *
 *          A constant (rigid) workspace disturbance D is supported via
 *          set_disturbance(): T_active(s) = D * T_ref(s). For a rigid D the body
 *          tangent tau(s) is invariant, so D only shifts the initial error e_0.
 */

#ifndef SERIAL_LINK_RMPCC_H
#define SERIAL_LINK_RMPCC_H

#include <Control/SerialLinkVelocityBase.h>
#include <Trajectory/CartesianSpline.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <memory>

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
    int    horizonSteps              = 12;       ///< Prediction horizon length N
    double terminalMultiplier        = 2.2;      ///< Extra weight on the last horizon stage
    bool   autoProgressRate          = true;     ///< Derive progressRateRef from trajectory duration
    double progressRateRef           = 0.0;      ///< Nominal progress rate sdot_ref (1/s)
    double autoProgressRateMin       = 0.0;      ///< Clamp for the auto-derived sdot_ref
    double autoProgressRateMax       = 1.15;     ///< Clamp for the auto-derived sdot_ref
    double progressRateMin           = 0.0;      ///< Lower bound on sdot
    double progressRateMax           = 0.0;      ///< Upper bound on sdot (0 => derive from ref)
    double progressRateMaxMultiplier = 1.0;      ///< progressRateMax = ref * multiplier when not set
    double progressUpperSlack        = 1e-4;     ///< Slack on the first-step terminal overrun constraint
    double progressScheduleSlack     = 0.0;      ///< Slack on the optional wall-clock schedule limit
    double completionTolerance       = 1e-3;     ///< Treat as "at the end" within this remaining progress

    // Numerics
    double tangentStep               = 1e-3;     ///< Finite-difference step for tau(s)
    double hessianRegularization     = 1e-8;     ///< Tikhonov term added to the QP Hessian

    // Progress shaping
    double progressReward            = 0.01;     ///< Linear reward encouraging forward progress
    double progressRateWeight        = 0.0;      ///< Quadratic tracking of sdot toward sdot_ref
    double progressRateSmoothWeight  = 0.5;      ///< Penalises sdot changes between stages

    // Riemannian tracking
    double lagWeight                 = 80.0;     ///< Weight on the along-path (lag) error component
    bool   poseFeedbackEnable        = true;     ///< Add an outer pose-feedback term after conversion to base-frame twist

    Eigen::Matrix<double, TWIST_DIM, TWIST_DIM> poseFeedbackGain =
        (Eigen::Vector<double, TWIST_DIM>() << 20.0, 20.0, 20.0, 5.0, 5.0, 5.0).finished().asDiagonal();
    Eigen::Matrix<double, TWIST_DIM, TWIST_DIM> metric =
        (Eigen::Vector<double, TWIST_DIM>() << 450.0, 450.0, 250.0, 10.0, 10.0, 10.0).finished().asDiagonal();
    Eigen::Matrix<double, TWIST_DIM, TWIST_DIM> contourWeight =
        (Eigen::Vector<double, TWIST_DIM>() << 600.0, 600.0, 350.0, 600.0, 600.0, 350.0).finished().asDiagonal();
    Eigen::Matrix<double, TWIST_DIM, TWIST_DIM> controlWeight =
        1e-5 * Eigen::Matrix<double, TWIST_DIM, TWIST_DIM>::Identity();
    Eigen::Matrix<double, TWIST_DIM, TWIST_DIM> controlRateWeight =
        1e-5 * Eigen::Matrix<double, TWIST_DIM, TWIST_DIM>::Identity();
    Eigen::Matrix<double, TWIST_DIM, TWIST_DIM> pathVelocityWeight =
        1e-5 * Eigen::Matrix<double, TWIST_DIM, TWIST_DIM>::Identity();

    Eigen::Vector3d linearVelocityMax  = Eigen::Vector3d(1.45, 0.34, 0.28);
    Eigen::Vector3d angularVelocityMax = Eigen::Vector3d(0.90, 0.90, 1.10);
};

/**
 * @brief Per-step diagnostics produced by a SerialLinkRMPCC control step.
 */
struct RmpccDiagnostics
{
    Eigen::Vector<double, 6> bodyTwist   = Eigen::Vector<double, 6>::Zero(); ///< Optimal first body twist u_1
    double progressRate          = 0.0;  ///< Optimal first progress rate sdot_1
    double pathProgress          = 0.0;  ///< Current externally estimated progress s
    double se3ErrorNorm          = 0.0;  ///< ||e_0|| (mixed units, diagnostic only)
    double contourError          = 0.0;  ///< Norm of the contour (perpendicular) error component
    double lagError              = 0.0;  ///< Signed lag (along-path) error component
    double rotationError         = 0.0;  ///< Rotation error magnitude (rad)
    double referenceLinearSpeed  = 0.0;  ///< |tau_lin| * sdot_ref
    double referenceAngularSpeed = 0.0;  ///< |tau_ang| * sdot_ref
    double qpStatus              = 1.0;  ///< 1 = solver succeeded, 0 = fallback used
    bool   fallbackUsed          = false;///< True when the QP solve failed and warm start was reused
};

/**
 * @brief Riemannian MPCC controller using a stored path and externally estimated progress.
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
         *         and an externally estimated progress.
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

        /**
         * @brief Set the current rigid workspace disturbance D (base frame). Default identity.
         *        T_active(s) = D * T_ref(s). Only shifts the initial error for rigid D.
         */
        void
        set_disturbance(const Eigen::Matrix4d &disturbance);

        /**
         * @brief Optional wall-clock schedule limit: the first step may not advance
         *        progress beyond this value (+ slack). Pass 1.0 to disable.
         */
        void
        set_schedule_limit(double scheduleProgressLimit);

        /**
         * @brief Reset progress, warm start and rate memory (call on a new goal).
         */
        void
        reset();

        /**
         * @brief Most recent externally supplied path progress s in [0, 1].
         */
        double
        path_progress() const { return _pathProgress; }

        /**
         * @brief Reference pose at a given progress, with the current disturbance applied.
         *        Non-const because CartesianSpline::query_state() mutates spline cache state.
         */
        RobotLibrary::Model::Pose
        reference_pose(double progress);

        /**
         * @brief Run one RMPCC control step from an externally estimated closest progress.
         * @param dt Time since the previous step [s].
         * @param estimatedProgress Current closest path progress in [0,1].
         */
        Eigen::VectorXd
        step(double dt, double estimatedProgress);

        /**
         * @brief Diagnostics from the most recent step().
         */
        const RmpccDiagnostics&
        diagnostics() const { return _diagnostics; }

    private:
        static constexpr int TWIST_DIM = 6;

        RmpccParameters _rmpcc;

        QPSolver<double> _qpSolver;

        RobotLibrary::Trajectory::CartesianSpline _trajectory;
        bool _trajectorySet = false;

        Eigen::Matrix4d _disturbance = Eigen::Matrix4d::Identity();
        double _scheduleProgressLimit = 1.0;

        // The external estimate is authoritative; predicted progress exists only in the QP.
        double _pathProgress = 0.0;
        double _lastProgressRate = 0.0;
        Eigen::Vector<double, TWIST_DIM> _lastBodyTwist = Eigen::Vector<double, TWIST_DIM>::Zero();
        Eigen::VectorXd _warmStart;

        RmpccDiagnostics _diagnostics;

        /**
         * @brief Reference transform at a progress value, with disturbance applied.
         *        Nominal path geometry comes from the trajectory; disturbance is applied here.
         *        Non-const: queries the (mutable) spline cache.
         */
        Eigen::Matrix4d
        reference_transform(double progress);

        /**
         * @brief Full reference body twist tau(s) * sdot_ref, using the trajectory's tangent.
         */
        Eigen::Vector<double, TWIST_DIM>
        body_twist_reference_at_progress(double progress);

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
         * @brief Build and solve the RMPCC QP; returns the optimal body twist and
         *        progress rate, and fills _diagnostics.
         * @param currentTransform Current endpoint pose as a 4x4 matrix.
         * @param dt               Control step [s].
         */
        void
        solve_rmpcc(const Eigen::Matrix4d &currentTransform, double dt);
};

} } // namespace

#endif
