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
#include <vector>

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

struct MpccObjectiveParameters
{
    double contourWeight = 1000.0;
    double lagWeight = 10.0;
    double orientationWeight = 40.0;
    double inputLinearWeight = 3e-4;
    double inputAngularWeight = 3e-4;
    double inputProgressWeight = 5e-4;
    double inputDifferenceWeight = 3e-5;
    double pathVelocityWeight = 2e-4;
    double progressReward = 2e-4;
};

struct MpccDiagnostics
{
    Eigen::Vector<double,6> error = Eigen::Vector<double,6>::Zero();
    Eigen::Vector<double,6> bodyTwist = Eigen::Vector<double,6>::Zero();
    Eigen::Vector<double,6> commandedBaseTwist = Eigen::Vector<double,6>::Zero();
    Eigen::Vector<double,6> realizedBaseTwist = Eigen::Vector<double,6>::Zero();
    Eigen::Vector<double,6> realizedBodyTwist = Eigen::Vector<double,6>::Zero();
    double progressRate = 0.0;
    double referenceProgress = 0.0;
    double nextProgress = 0.0;
    double positionError = 0.0;
    double orientationError = 0.0;
    double qpStatus = 1.0;
    bool qpConverged = false;
    double qpIterations = 0.0;
    double qpFinalStepSize = 0.0;
    double qpMaximumConstraintViolation = 0.0;
    double qpObjective = 0.0;
    double qpSolveTimeSeconds = 0.0;
    double twistRealizationError = 0.0;
    bool parentFrameMotionActive = false;
    Eigen::Vector<double,6> parentFrameBodyTwist = Eigen::Vector<double,6>::Zero();
    Eigen::Matrix4d measuredParentPose = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d predictedParentPoseFirst = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d predictedParentPoseHorizon = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d parentReferenceFactorFirst = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d repairedReferenceDisplacementFirst = Eigen::Matrix4d::Identity();
    std::vector<Eigen::Matrix4d> predictedParentTransforms;
    double parentMeasurementTimeSeconds = 0.0;
    double parentMeasurementAgeSeconds = 0.0;
    int stageControlDimension = 7;
    Eigen::VectorXd optimalHorizon;
    /** Number of independent equality rows carried by the frozen QP. */
    double qpEqualityRows = 0.0;
    /** Maximum |A x - y| over those equality rows for the returned solution. */
    double qpEqualityViolation = 0.0;
    /** Exact frozen condensed QP used by the most recent solve (diagnostic only). */
    Eigen::MatrixXd qpHessian;
    Eigen::VectorXd qpGradient;
    Eigen::MatrixXd qpConstraintMatrix;
    Eigen::VectorXd qpConstraintVector;
    /** Fixed-variable equality block. Empty when no bound pair coincides. */
    Eigen::MatrixXd qpEqualityMatrix;
    Eigen::VectorXd qpEqualityVector;
    Eigen::VectorXd qpSeed;
    /** Last vector returned by the QP solver, including a failing solve. */
    Eigen::VectorXd qpReturnedSolution;
    Eigen::VectorXd shiftedWarmStart;
};

/** Read-only horizon geometry passed to optional MPCC QP extensions. */
struct MpccQpExtensionContext
{
    int horizon = 0;
    int baseControlDimension = 0;
    int stageControlDimension = 7;
    double dt = 0.0;
    Eigen::Matrix3d predictionRotation = Eigen::Matrix3d::Identity();
    Eigen::Vector3d currentToolPositionBase = Eigen::Vector3d::Zero();
    Eigen::Matrix3d currentToolRotationBase = Eigen::Matrix3d::Identity();
    std::vector<Eigen::Matrix4d> parentTransforms;
    std::vector<double> stageProgress;
    Eigen::VectorXd previousWarmStart;
    /**
     * Whether the QP reached its convergence exit rather than the iteration
     * cap. Meaningful only in on_extended_qp_solution(); it is false while the
     * problem is still being assembled in extend_qp_problem().
     */
    bool qpConverged = false;
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

        /** Configure the condensed-QP objective without changing its structure. */
        void
        set_objective_parameters(const MpccObjectiveParameters &parameters);

        /** Explicit normalized progress-rate bounds, callable after set_trajectory(). */
        void
        set_progress_rate_limits(double minimum, double nominal, double maximum);

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

        /** Run a step while compensating a timestamped parent-frame sample to control time. */
        Eigen::VectorXd
        step_at_time(double controlTimeSeconds, double dt);

        /** Timestamp-compensated step with an externally supplied initial progress. */
        Eigen::VectorXd
        step_at_time(double controlTimeSeconds, double dt, double estimatedProgress);

        /**
         * @brief Current internally integrated path progress.
         */
        double
        path_progress() const { return _pathProgress; }

        const MpccDiagnostics &
        diagnostics() const { return _diagnostics; }

    protected:
        static constexpr int ERROR_DIM = 6;
        static constexpr int PHYSICAL_TWIST_DIMENSION = 6;
        static constexpr int PROGRESS_RATE_INDEX = 6;
        static constexpr int BASE_STAGE_CONTROL_DIMENSION = 7;

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

        const RobotLibrary::Trajectory::CartesianSpline &
        reference_trajectory() const { return _trajectory; }

        const RobotLibrary::Trajectory::CartesianTrajectoryFrameState &
        trajectory_frame() const { return _trajectoryFrame; }

        /** Reference pose used for the current tracking error. */
        virtual RobotLibrary::Model::Pose
        reference_pose_at_progress(double progress);

        /** Number of decision variables in each horizon stage. */
        virtual int
        stage_control_dimension() const { return BASE_STAGE_CONTROL_DIMENSION; }

        /**
         * Reference-coordinate tangents for virtual inputs after s_dot.
         * Columns are expressed in the same fixed local prediction frame as
         * the physical twist and path tangent.  The predictor inserts them as
         * -dt*tangent*virtual_rate.
         */
        virtual Eigen::MatrixXd
        additional_reference_tangents(
            int stage,
            double progress,
            const Eigen::Matrix3d &predictionRotation,
            const Eigen::Matrix4d &predictedParentTransform) const;

        /** Configure box bounds and a feasible nominal value for extra inputs. */
        virtual void
        configure_additional_stage_inputs(int stage,
                                          Eigen::VectorXd &lower,
                                          Eigen::VectorXd &upper,
                                          Eigen::VectorXd &nominal) const;

        /** Reset controller-owned virtual coordinates on a new trajectory. */
        virtual void
        reset_additional_virtual_state();

        /** Default returns the original MPCC position weight unchanged. */
        virtual Eigen::Matrix3d
        position_error_weight(
            int stage,
            const Eigen::Vector<double,ERROR_DIM> &pathTangent,
            const Eigen::Matrix3d &defaultWeight,
            const Eigen::Matrix3d &predictionRotation,
            const Eigen::Matrix4d &predictedParentTransform) const;

        /** Optional cost/constraint/variable augmentation; default is a no-op. */
        virtual void
        extend_qp_problem(const MpccQpExtensionContext &context,
                          Eigen::MatrixXd &hessian,
                          Eigen::VectorXd &gradient,
                          Eigen::MatrixXd &constraintMatrix,
                          Eigen::VectorXd &constraintVector,
                          Eigen::VectorXd &seed);

        /** Shift variables appended after the stage-stacked decision vector. */
        virtual void
        shift_extension_warm_start(const MpccQpExtensionContext &context,
                                   const Eigen::VectorXd &optimum,
                                   Eigen::VectorXd &shiftedWarmStart);

        /** Observe a valid augmented-QP solution. */
        virtual void
        on_extended_qp_solution(const MpccQpExtensionContext &context,
                                const Eigen::VectorXd &optimum);

        /** Optional composition point between MPCC and resolved-rate IK. */
        virtual Eigen::Vector<double,6>
        postprocess_base_twist(const Eigen::Vector<double,6> &baseTwist,
                               double dt);

        /** Observe the endpoint twist actually realizable by the joint command. */
        virtual void
        on_twist_resolved(const Eigen::Vector<double,6> &commandedBaseTwist,
                          const Eigen::Vector<double,6> &realizedBaseTwist,
                          double dt);

    private:
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
        Eigen::Vector<double,BASE_STAGE_CONTROL_DIMENSION> _uLast =
            Eigen::Vector<double,BASE_STAGE_CONTROL_DIMENSION>::Zero();
        Eigen::VectorXd _warmStart;

        QPSolver<double> _qpSolver;
        double _qpStepSizeTolerance = 1e-5;
        MpccDiagnostics _diagnostics;

        Eigen::Vector<double,BASE_STAGE_CONTROL_DIMENSION>
        solve_mpcc(const Eigen::Vector<double,ERROR_DIM> &error0,
                   const Eigen::Matrix3d &referenceRotation,
                   const Eigen::Vector3d &currentToolPositionBase,
                   const Eigen::Matrix3d &currentToolRotationBase,
                   double parentMeasurementAgeSeconds);

        Eigen::VectorXd
        step_impl(double dt, double parentMeasurementAgeSeconds);
};

} } // namespace

#endif
