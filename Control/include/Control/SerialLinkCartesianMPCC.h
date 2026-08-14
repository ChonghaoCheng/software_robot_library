/**
 * @file    SerialLinkCartesianMPCC.h
 * @brief   Cartesian position MPCC with independent rotation-error tracking.
 */

#ifndef SERIAL_LINK_CARTESIAN_MPCC_H
#define SERIAL_LINK_CARTESIAN_MPCC_H

#include <Control/SerialLinkMPCC.h>

namespace RobotLibrary { namespace Control {

/**
 * @brief MPCC variant using a Cartesian position tangent and a separate
 *        rotation-error channel.
 *
 * The position tangent is formed from the Euclidean Cartesian path derivative.
 * Orientation remains a separate rotation-error channel, but its reference
 * angular derivative is retained in the progress-coupled prediction so a
 * rotating path has the required angular feedforward.
 *
 * This class and SerialLinkMPCC currently share the same local linear error
 * predictor and QP. Direct Lie-group horizon propagation is implemented by
 * SerialLinkRMPCC, not by SerialLinkMPCC.
 */
class SerialLinkCartesianMPCC : public SerialLinkMPCC
{
    public:
        using SerialLinkMPCC::SerialLinkMPCC;

    protected:
        /**
         * @brief Cartesian position tangent plus the independent SO(3)
         *        reference derivative, both expressed in the reference pose.
         *
         * Exposed as a protected static helper so its rotating-reference
         * semantics can be tested without constructing a robot model.
         */
        static Eigen::Vector<double,6>
        cartesian_path_tangent(
            RobotLibrary::Trajectory::CartesianSpline &trajectory,
            double progress);

        Eigen::Vector<double,6>
        path_tangent_at_progress(double progress,
                                 const Eigen::Matrix3d &referenceRotation) override;
};

} } // namespace RobotLibrary::Control

#endif
