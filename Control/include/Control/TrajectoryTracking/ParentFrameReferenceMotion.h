/**
 * @file ParentFrameReferenceMotion.h
 * @brief Causal parent-frame prediction and stage reference-motion factors.
 */

#ifndef CONTROL_TRAJECTORY_TRACKING_PARENT_FRAME_REFERENCE_MOTION_H
#define CONTROL_TRAJECTORY_TRACKING_PARENT_FRAME_REFERENCE_MOTION_H

#include <Math/MathFunctions.h>
#include <Trajectory/CartesianTrajectoryFrame.h>

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace RobotLibrary { namespace Control {

/**
 * Causal constant-body-twist model of a trajectory parent frame A(t).
 * The estimate uses only the two latest measured poses and their timestamps.
 */
class CausalParentFrameMotion
{
    public:
        enum class UpdateStatus
        {
            NeverUpdated = 0,
            StaticPose = 1,
            FirstTimestampedSample = 2,
            TimestampAdvanced = 3,
            DuplicateTimestampIgnored = 4,
            OutOfOrderTimestampIgnored = 5
        };

        void reset()
        {
            _current = Eigen::Matrix4d::Identity();
            _previous = Eigen::Matrix4d::Identity();
            _bodyTwist.setZero();
            _currentTime = 0.0;
            _previousTime = 0.0;
            _hasCurrent = false;
            _hasPrevious = false;
            _hasVelocity = false;
            _lastStatus = UpdateStatus::NeverUpdated;
            _timestampedSetterCount = 0;
            _staticSetterCount = 0;
            _duplicateTimestampCount = 0;
            _outOfOrderTimestampCount = 0;
            _lastElapsed = 0.0;
            _minimumPositiveElapsed = std::numeric_limits<double>::infinity();
            _maximumPositiveElapsed = 0.0;
        }

        void set_static_pose(const Eigen::Matrix4d &pose)
        {
            require_finite_pose(pose);
            _current = pose;
            _previous = pose;
            _bodyTwist.setZero();
            _currentTime = 0.0;
            _previousTime = 0.0;
            _hasCurrent = true;
            _hasPrevious = false;
            _hasVelocity = false;
            ++_staticSetterCount;
            _lastStatus = UpdateStatus::StaticPose;
        }

        void update(const Eigen::Matrix4d &pose, const double timestampSeconds)
        {
            ++_timestampedSetterCount;
            require_finite_pose(pose);
            if(not std::isfinite(timestampSeconds))
            {
                throw std::invalid_argument(
                    "[ERROR] [PARENT FRAME MOTION] Timestamp must be finite.");
            }

            if(_hasCurrent and timestampSeconds > _currentTime)
            {
                _previous = _current;
                _previousTime = _currentTime;
                _hasPrevious = true;
                _current = pose;
                _currentTime = timestampSeconds;
                const double elapsed = _currentTime - _previousTime;
                _lastElapsed = elapsed;
                _minimumPositiveElapsed = std::min(_minimumPositiveElapsed, elapsed);
                _maximumPositiveElapsed = std::max(_maximumPositiveElapsed, elapsed);
                _bodyTwist = RobotLibrary::Math::se3_logarithm(
                    RobotLibrary::Math::se3_inverse(_previous) * _current) / elapsed;
                _hasVelocity = true;
                _lastStatus = UpdateStatus::TimestampAdvanced;
                return;
            }

            // A repeated/out-of-order measurement is not a new sample. Keep
            // the latest two-sample estimate rather than creating a zero/spike
            // sequence when TF is published below the controller frequency.
            if(_hasCurrent)
            {
                if(timestampSeconds == _currentTime)
                {
                    ++_duplicateTimestampCount;
                    _lastStatus = UpdateStatus::DuplicateTimestampIgnored;
                }
                else
                {
                    ++_outOfOrderTimestampCount;
                    _lastStatus = UpdateStatus::OutOfOrderTimestampIgnored;
                }
                return;
            }

            // The first sample establishes A_k but cannot define velocity.
            _current = pose;
            _currentTime = timestampSeconds;
            _bodyTwist.setZero();
            _hasVelocity = false;
            _hasCurrent = true;
            _hasPrevious = false;
            _lastStatus = UpdateStatus::FirstTimestampedSample;
        }

        bool has_velocity() const { return _hasVelocity; }

        const Eigen::Vector<double,6> &body_twist() const { return _bodyTwist; }

        const Eigen::Matrix4d &current_pose() const { return _current; }

        double current_time() const { return _currentTime; }
        double previous_time() const { return _previousTime; }
        double last_elapsed() const { return _lastElapsed; }
        double minimum_positive_elapsed() const
        {
            return std::isfinite(_minimumPositiveElapsed) ? _minimumPositiveElapsed : 0.0;
        }
        double maximum_positive_elapsed() const { return _maximumPositiveElapsed; }
        UpdateStatus last_update_status() const { return _lastStatus; }
        std::uint64_t timestamped_setter_count() const { return _timestampedSetterCount; }
        std::uint64_t static_setter_count() const { return _staticSetterCount; }
        std::uint64_t duplicate_timestamp_count() const { return _duplicateTimestampCount; }
        std::uint64_t out_of_order_timestamp_count() const { return _outOfOrderTimestampCount; }
        bool too_small_interval_observable_without_policy_change() const { return false; }

        Eigen::Matrix4d predicted_pose(const int stage, const double dt) const
        {
            if(stage < 0 or not std::isfinite(dt) or dt <= 0.0)
            {
                throw std::invalid_argument(
                    "[ERROR] [PARENT FRAME MOTION] Nonnegative stage and positive dt required.");
            }
            if(not _hasCurrent)
            {
                return Eigen::Matrix4d::Identity();
            }
            return _current * RobotLibrary::Math::se3_exponential(
                static_cast<double>(stage) * dt * _bodyTwist);
        }

    private:
        static void require_finite_pose(const Eigen::Matrix4d &pose)
        {
            if(not pose.allFinite())
            {
                throw std::invalid_argument(
                    "[ERROR] [PARENT FRAME MOTION] Parent pose must be finite.");
            }
        }

        Eigen::Matrix4d _current = Eigen::Matrix4d::Identity();
        Eigen::Matrix4d _previous = Eigen::Matrix4d::Identity();
        Eigen::Vector<double,6> _bodyTwist = Eigen::Vector<double,6>::Zero();
        double _currentTime = 0.0;
        double _previousTime = 0.0;
        bool _hasCurrent = false;
        bool _hasPrevious = false;
        bool _hasVelocity = false;
        UpdateStatus _lastStatus = UpdateStatus::NeverUpdated;
        std::uint64_t _timestampedSetterCount = 0;
        std::uint64_t _staticSetterCount = 0;
        std::uint64_t _duplicateTimestampCount = 0;
        std::uint64_t _outOfOrderTimestampCount = 0;
        double _lastElapsed = 0.0;
        double _minimumPositiveElapsed = std::numeric_limits<double>::infinity();
        double _maximumPositiveElapsed = 0.0;
};

/**
 * Predicted parent-frame state for a horizon stage.
 *
 * CausalParentFrameMotion stores the parent twist in parent/body axes.  The
 * trajectory-frame contract requires the frame-origin point twist in base
 * axes, so both linear and angular components are rotated by the predicted
 * stage orientation before moving-frame transport is applied.
 */
inline RobotLibrary::Trajectory::CartesianTrajectoryFrameState
predicted_parent_frame_state(const CausalParentFrameMotion &motion,
                             const int stage,
                             const double dt)
{
    RobotLibrary::Trajectory::CartesianTrajectoryFrameState state;
    state.transformInBase = motion.predicted_pose(stage, dt);
    const Eigen::Matrix3d rotation = state.transformInBase.block<3,3>(0,0);
    state.twistInBase.head<3>() = rotation * motion.body_twist().head<3>();
    state.twistInBase.tail<3>() = rotation * motion.body_twist().tail<3>();
    state.measurementTimeSeconds = motion.current_time();
    return state;
}

inline Eigen::Matrix4d
parent_frame_reference_factor(const Eigen::Matrix4d &pathPose,
                              const Eigen::Matrix4d &parentCurrent,
                              const Eigen::Matrix4d &parentNext)
{
    return RobotLibrary::Math::se3_inverse(pathPose)
           * RobotLibrary::Math::se3_inverse(parentCurrent)
           * parentNext * pathPose;
}

inline Eigen::Matrix4d
legacy_repaired_reference_displacement(
    const Eigen::Matrix4d &pathPose,
    const Eigen::Vector<double,6> &pathTangent,
    const double progressRate,
    const double dt,
    const Eigen::Matrix4d &parentCurrent,
    const Eigen::Matrix4d &parentNext)
{
    return parent_frame_reference_factor(pathPose, parentCurrent, parentNext)
           * RobotLibrary::Math::se3_exponential(dt * pathTangent * progressRate);
}

} } // namespace RobotLibrary::Control

#endif
