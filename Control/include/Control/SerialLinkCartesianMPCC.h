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
 * Unlike SerialLinkMPCC, this controller does not use the rotational part of
 * an SE(3) path tangent in the progress-coupled prediction. The position
 * tangent drives contouring/lag control; orientation is handled by the
 * existing rotation-error term.
 */
class SerialLinkCartesianMPCC : public SerialLinkMPCC
{
    public:
        using SerialLinkMPCC::SerialLinkMPCC;

    protected:
        Eigen::Vector<double,6>
        path_tangent_at_progress(double progress,
                                 const Eigen::Matrix3d &referenceRotation) override;
};

} } // namespace RobotLibrary::Control

#endif
