/**
 * @file SerialLinkTimeIndexedMPC.h
 * @brief Shared interface for time-indexed Cartesian MPC controllers.
 */

#ifndef SERIAL_LINK_TIME_INDEXED_MPC_H
#define SERIAL_LINK_TIME_INDEXED_MPC_H

#include <Control/Core/SerialLinkVelocityBase.h>
#include <Trajectory/CartesianSpline.h>
#include <Trajectory/CartesianTrajectoryFrame.h>

namespace RobotLibrary { namespace Control {

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
};

} } // namespace RobotLibrary::Control

#endif
