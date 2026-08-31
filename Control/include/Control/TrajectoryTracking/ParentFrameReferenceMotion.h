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
            NoNewMeasurement = 4,
            DuplicateTimestampPoseMismatchRejected = 5,
            OutOfOrderTimestampRejected = 6,
            TooSmallIntervalReseeded = 7,
            GenerationResetFirstSample = 8,
            NonfiniteInputRejected = 9
        };

        static constexpr double minimumElapsedSeconds = 1e-6;
        static constexpr double maximumMeasurementAgeSeconds = 0.050;
        static constexpr double duplicateTranslationTolerance = 1e-10;
        static constexpr double duplicateRotationTolerance = 1e-10;

        void reset()
        {
            _current = Eigen::Matrix4d::Identity();
            _previous = Eigen::Matrix4d::Identity();
            _bodyTwist.setZero();
            _currentTime = 0.0;
            _previousTime = 0.0;
            _hasCurrent = false;
            _hasPrevious = false;
            _hasTimestampedCurrent = false;
            _rawVelocityValid = false;
            _hasGeneration = false;
            _generation = 0;
            _evaluationTime = 0.0;
            _measurementAge = 0.0;
            _lastStatus = UpdateStatus::NeverUpdated;
            _timestampedSetterCount = 0;
            _staticSetterCount = 0;
            _duplicateTimestampCount = 0;
            _duplicatePoseMismatchCount = 0;
            _outOfOrderTimestampCount = 0;
            _tooSmallIntervalCount = 0;
            _generationResetCount = 0;
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
            _hasTimestampedCurrent = false;
            _rawVelocityValid = false;
            _hasGeneration = false;
            _measurementAge = 0.0;
            ++_staticSetterCount;
            _lastStatus = UpdateStatus::StaticPose;
        }

        UpdateStatus update(const Eigen::Matrix4d &pose,
                            const double timestampSeconds)
        {
            return update_impl(pose, timestampSeconds, 0, timestampSeconds, false);
        }

        UpdateStatus update(const Eigen::Matrix4d &pose,
                            const double timestampSeconds,
                            const std::uint64_t generation,
                            const double evaluationTimeSeconds)
        {
            return update_impl(
                pose, timestampSeconds, generation, evaluationTimeSeconds, true);
        }

        void evaluate_at(const double evaluationTimeSeconds)
        {
            if(not std::isfinite(evaluationTimeSeconds))
            {
                _lastStatus = UpdateStatus::NonfiniteInputRejected;
                throw std::invalid_argument(
                    "[ERROR] [PARENT FRAME MOTION] Evaluation time must be finite.");
            }
            _evaluationTime = evaluationTimeSeconds;
            update_measurement_age();
        }

        bool has_velocity() const { return _rawVelocityValid; }
        bool raw_velocity_valid() const { return _rawVelocityValid; }
        bool prediction_velocity_active() const
        {
            return _rawVelocityValid
                && _measurementAge <= maximumMeasurementAgeSeconds;
        }

        const Eigen::Vector<double,6> &body_twist() const { return _bodyTwist; }

        Eigen::Vector<double,6> prediction_body_twist() const
        {
            return prediction_velocity_active()
                ? _bodyTwist : Eigen::Vector<double,6>::Zero();
        }

        const Eigen::Matrix4d &current_pose() const { return _current; }

        double current_time() const { return _currentTime; }
        double previous_time() const { return _previousTime; }
        double evaluation_time() const { return _evaluationTime; }
        double measurement_age() const { return _measurementAge; }
        double last_elapsed() const { return _lastElapsed; }
        double minimum_positive_elapsed() const
        {
            return std::isfinite(_minimumPositiveElapsed) ? _minimumPositiveElapsed : 0.0;
        }
        double maximum_positive_elapsed() const { return _maximumPositiveElapsed; }
        UpdateStatus last_update_status() const { return _lastStatus; }
        std::uint64_t current_generation() const { return _generation; }
        bool has_generation() const { return _hasGeneration; }
        std::uint64_t timestamped_setter_count() const { return _timestampedSetterCount; }
        std::uint64_t static_setter_count() const { return _staticSetterCount; }
        std::uint64_t duplicate_timestamp_count() const { return _duplicateTimestampCount; }
        std::uint64_t duplicate_pose_mismatch_count() const
        {
            return _duplicatePoseMismatchCount;
        }
        std::uint64_t out_of_order_timestamp_count() const { return _outOfOrderTimestampCount; }
        std::uint64_t too_small_interval_count() const { return _tooSmallIntervalCount; }
        std::uint64_t generation_reset_count() const { return _generationResetCount; }
        bool too_small_interval_observable_without_policy_change() const { return true; }

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
                static_cast<double>(stage) * dt * prediction_body_twist());
        }

    private:
        UpdateStatus update_impl(const Eigen::Matrix4d &pose,
                                 const double timestampSeconds,
                                 const std::uint64_t generation,
                                 const double evaluationTimeSeconds,
                                 const bool generationProvided)
        {
            ++_timestampedSetterCount;
            if(not pose.allFinite() || not std::isfinite(timestampSeconds)
               || not std::isfinite(evaluationTimeSeconds))
            {
                _lastStatus = UpdateStatus::NonfiniteInputRejected;
                throw std::invalid_argument(
                    "[ERROR] [PARENT FRAME MOTION] Pose and timestamps must be finite.");
            }
            _evaluationTime = evaluationTimeSeconds;

            if(generationProvided && _hasGeneration && generation != _generation)
            {
                ++_generationResetCount;
                _generation = generation;
                seed_timestamped_sample(pose, timestampSeconds);
                _lastStatus = UpdateStatus::GenerationResetFirstSample;
                update_measurement_age();
                return _lastStatus;
            }
            if(generationProvided && !_hasGeneration)
            {
                _generation = generation;
                _hasGeneration = true;
            }

            if(!_hasTimestampedCurrent)
            {
                seed_timestamped_sample(pose, timestampSeconds);
                _lastStatus = UpdateStatus::FirstTimestampedSample;
                update_measurement_age();
                return _lastStatus;
            }

            if(timestampSeconds == _currentTime)
            {
                if(same_pose_within_duplicate_tolerance(pose, _current))
                {
                    ++_duplicateTimestampCount;
                    _lastStatus = UpdateStatus::NoNewMeasurement;
                }
                else
                {
                    ++_duplicatePoseMismatchCount;
                    _lastStatus = UpdateStatus::DuplicateTimestampPoseMismatchRejected;
                }
                update_measurement_age();
                return _lastStatus;
            }

            if(timestampSeconds < _currentTime)
            {
                ++_outOfOrderTimestampCount;
                _lastStatus = UpdateStatus::OutOfOrderTimestampRejected;
                update_measurement_age();
                return _lastStatus;
            }

            const double elapsed = timestampSeconds - _currentTime;
            if(elapsed < minimumElapsedSeconds)
            {
                ++_tooSmallIntervalCount;
                seed_timestamped_sample(pose, timestampSeconds);
                _lastStatus = UpdateStatus::TooSmallIntervalReseeded;
                update_measurement_age();
                return _lastStatus;
            }

            _previous = _current;
            _previousTime = _currentTime;
            _hasPrevious = true;
            _current = pose;
            _currentTime = timestampSeconds;
            _lastElapsed = elapsed;
            _minimumPositiveElapsed = std::min(_minimumPositiveElapsed, elapsed);
            _maximumPositiveElapsed = std::max(_maximumPositiveElapsed, elapsed);
            _bodyTwist = RobotLibrary::Math::se3_logarithm(
                RobotLibrary::Math::se3_inverse(_previous) * _current) / elapsed;
            _rawVelocityValid = true;
            _lastStatus = UpdateStatus::TimestampAdvanced;
            update_measurement_age();
            return _lastStatus;
        }

        void seed_timestamped_sample(const Eigen::Matrix4d &pose,
                                     const double timestampSeconds)
        {
            _current = pose;
            _previous = pose;
            _currentTime = timestampSeconds;
            _previousTime = timestampSeconds;
            _bodyTwist.setZero();
            _hasCurrent = true;
            _hasPrevious = false;
            _hasTimestampedCurrent = true;
            _rawVelocityValid = false;
            _lastElapsed = 0.0;
        }

        void update_measurement_age()
        {
            _measurementAge = _hasTimestampedCurrent
                ? std::max(0.0, _evaluationTime - _currentTime) : 0.0;
        }

        static bool same_pose_within_duplicate_tolerance(
            const Eigen::Matrix4d &a, const Eigen::Matrix4d &b)
        {
            const double translation =
                (a.block<3,1>(0,3) - b.block<3,1>(0,3)).norm();
            const Eigen::Matrix3d relative =
                b.block<3,3>(0,0).transpose() * a.block<3,3>(0,0);
            const double cosine = std::clamp((relative.trace() - 1.0) * 0.5, -1.0, 1.0);
            return translation <= duplicateTranslationTolerance
                && std::acos(cosine) <= duplicateRotationTolerance;
        }

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
        bool _hasTimestampedCurrent = false;
        bool _rawVelocityValid = false;
        bool _hasGeneration = false;
        std::uint64_t _generation = 0;
        double _evaluationTime = 0.0;
        double _measurementAge = 0.0;
        UpdateStatus _lastStatus = UpdateStatus::NeverUpdated;
        std::uint64_t _timestampedSetterCount = 0;
        std::uint64_t _staticSetterCount = 0;
        std::uint64_t _duplicateTimestampCount = 0;
        std::uint64_t _duplicatePoseMismatchCount = 0;
        std::uint64_t _outOfOrderTimestampCount = 0;
        std::uint64_t _tooSmallIntervalCount = 0;
        std::uint64_t _generationResetCount = 0;
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
    const Eigen::Vector<double,6> predictionTwist = motion.prediction_body_twist();
    state.twistInBase.head<3>() = rotation * predictionTwist.head<3>();
    state.twistInBase.tail<3>() = rotation * predictionTwist.tail<3>();
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
