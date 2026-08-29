/**
 * @file SerialLinkTimeIndexedMPC.h
 * @brief Shared interface for time-indexed Cartesian MPC controllers.
 */

#ifndef SERIAL_LINK_TIME_INDEXED_MPC_H
#define SERIAL_LINK_TIME_INDEXED_MPC_H

#include <Control/Core/SerialLinkVelocityBase.h>
#include <Trajectory/CartesianSpline.h>
#include <Trajectory/CartesianTrajectoryFrame.h>

#include <Eigen/Core>

namespace RobotLibrary { namespace Control {

/**
 * @brief Read-only observability record for one time-indexed MPC invocation.
 *
 * Every field mirrors a quantity the controller already computed; recording
 * them changes no control decision. `qpStatus` is 1.0 on any invocation that
 * returns normally: the underlying solvers report failure by throwing, so a
 * returned invocation is by construction a solved one. Failures continue to
 * propagate through that exception path and are never reported here as
 * successes.
 */
struct TimeIndexedMpcDiagnostics
{
    double invocationIndex = 0.0;               ///< Monotonic count of completed invocations.
    double controllerDt = 0.0;                  ///< Prediction/integration dt used this invocation.
    double referenceTime = 0.0;                 ///< Clamped trajectory time the reference was queried at.
    double qpStatus = 1.0;                      ///< 1.0 = solved and returned; failure throws instead.
    double qpIterations = 0.0;
    double qpFinalStepSize = 0.0;
    double qpPrimalViolation = 0.0;
    bool qpConverged = false;
    bool qpHitMaxIterations = false;
    Eigen::Vector<double,6> commandedTwist =
        Eigen::Vector<double,6>::Zero();        ///< Post-clamp endpoint twist in base axes.
    Eigen::Vector<double,6> clampFrameTwist =
        Eigen::Vector<double,6>::Zero();        ///< Same twist in the axes the clamp was applied in.
    Eigen::Vector<double,6> realizedTwist =
        Eigen::Vector<double,6>::Zero();        ///< J * qdot from the same Jacobian and returned command.
    double twistRealizationError = 0.0;         ///< ||realizedTwist - commandedTwist||.
    bool linearLimitActive = false;             ///< A linear axis sits on the Cartesian clamp.
    bool angularLimitActive = false;            ///< An angular axis sits on the Cartesian clamp.
    bool jointLimitActive = false;              ///< A joint command sits on its control limit.
    double jointLimitMargin = 0.0;              ///< Smallest distance from any joint command to
                                                ///< its control limit; jointLimitActive is exactly
                                                ///< this margin falling to zero.
    double effectiveLinearLimit = 0.0;          ///< Per-axis linear clamp in force this invocation.
    double effectiveAngularLimit = 0.0;         ///< Per-axis angular clamp in force this invocation.
    double solveTimeSeconds = 0.0;              ///< Wall time spent inside the controller call.
    bool finite = true;                         ///< False if any recorded quantity is NaN/Inf.
};

/**
 * @brief Common action-server contract for time-indexed Cartesian MPCs.
 *
 * This interface deliberately contains no prediction model. Cartesian MPC and
 * Lie-algebra MPC share trajectory/frame plumbing and the resolved-rate layer,
 * while retaining independent state and dynamics definitions.
 */
class SerialLinkTimeIndexedMPC : public SerialLinkVelocityBase
{
    public:
        using SerialLinkVelocityBase::SerialLinkVelocityBase;
        virtual ~SerialLinkTimeIndexedMPC() = default;

        virtual void
        set_trajectory(const RobotLibrary::Trajectory::CartesianSpline &trajectory) = 0;

        virtual void
        set_trajectory_frame(
            const RobotLibrary::Trajectory::CartesianTrajectoryFrameState &frame) = 0;

        virtual void clear_trajectory() = 0;
        virtual bool has_trajectory() const = 0;

        virtual Eigen::VectorXd
        track_endpoint_trajectory_at_time(const double &time) = 0;

        /** @brief Diagnostics from the most recent successful invocation. */
        const TimeIndexedMpcDiagnostics &diagnostics() const { return _timeIndexedDiagnostics; }

    protected:
        /**
         * @brief Record one invocation's diagnostics. Behaviour-neutral.
         *
         * Computes the realized twist from the same Jacobian and joint command
         * the caller is about to return, so the recorded realization error is
         * the one actually incurred rather than a reconstruction.
         *
         * @param commandedTwistInBase Post-clamp endpoint twist handed to the resolved-rate layer.
         * @param clampFrameTwist      The same twist expressed in the axes the clamp was applied
         *                             in. Cartesian MPC clamps in base axes and passes the same
         *                             vector twice; Lie-algebra MPC clamps in body axes, where a
         *                             per-axis comparison against the limit is the meaningful one.
         * @param jointCommand         Joint velocities returned by resolve_endpoint_twist().
         * @param referenceTime        Clamped trajectory time used for the reference query.
         * @param controllerDt         Prediction/integration dt used this invocation.
         * @param maxLinearSpeed       Per-axis linear clamp applied to the twist.
         * @param maxAngularSpeed      Per-axis angular clamp applied to the twist.
         * @param solveTimeSeconds     Wall time spent in the controller call.
         */
        void
        record_time_indexed_diagnostics(
            const Eigen::Vector<double,6> &commandedTwistInBase,
            const Eigen::Vector<double,6> &clampFrameTwist,
            const Eigen::VectorXd &jointCommand,
            const double referenceTime,
            const double controllerDt,
            const double maxLinearSpeed,
            const double maxAngularSpeed,
            const double solveTimeSeconds);

        TimeIndexedMpcDiagnostics _timeIndexedDiagnostics;
};

} } // namespace RobotLibrary::Control

#endif
