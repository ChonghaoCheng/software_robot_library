/**
 * @file ParentFrameReferenceMotion.h
 * @brief Causal parent-frame prediction and stage reference-motion factors.
 */

#ifndef CONTROL_TRAJECTORY_TRACKING_PARENT_FRAME_REFERENCE_MOTION_H
#define CONTROL_TRAJECTORY_TRACKING_PARENT_FRAME_REFERENCE_MOTION_H

#include <Math/MathFunctions.h>

#include <Eigen/Core>

#include <cmath>
#include <stdexcept>

namespace RobotLibrary { namespace Control {

/**
 * Causal constant-body-twist model of a trajectory parent frame A(t).
 * The estimate uses only the two latest measured poses and their timestamps.
 */
class CausalParentFrameMotion
{
    public:
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
        }

        void update(const Eigen::Matrix4d &pose, const double timestampSeconds)
        {
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
                _bodyTwist = RobotLibrary::Math::se3_logarithm(
                    RobotLibrary::Math::se3_inverse(_previous) * _current) / elapsed;
                _hasVelocity = true;
                return;
            }

            // A repeated/out-of-order measurement is not a new sample. Keep
            // the latest two-sample estimate rather than creating a zero/spike
            // sequence when TF is published below the controller frequency.
            if(_hasCurrent)
            {
                return;
            }

            // The first sample establishes A_k but cannot define velocity.
            _current = pose;
            _currentTime = timestampSeconds;
            _bodyTwist.setZero();
            _hasVelocity = false;
            _hasCurrent = true;
            _hasPrevious = false;
        }

        bool has_velocity() const { return _hasVelocity; }

        const Eigen::Vector<double,6> &body_twist() const { return _bodyTwist; }

        const Eigen::Matrix4d &current_pose() const { return _current; }

        double current_time() const { return _currentTime; }

        /** Diagnostic-branch-only frozen-state counterfactual hook. */
        void diagnostic_override_body_twist(
            const Eigen::Vector<double,6> &bodyTwist)
        {
            if(not bodyTwist.allFinite())
                throw std::invalid_argument("Diagnostic parent twist must be finite.");
            _bodyTwist = bodyTwist;
            _hasVelocity = bodyTwist.squaredNorm() > 0.0;
        }

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
};

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
