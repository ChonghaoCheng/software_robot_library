/** @file SerialLinkAdmittanceMPCC.cpp */

#include <Control/Contact/SerialLinkAdmittanceMPCC.h>

namespace RobotLibrary { namespace Control {

Eigen::Vector<double,6>
SerialLinkAdmittanceMPCC::postprocess_base_twist(
    const Eigen::Vector<double,6> &baseTwist,
    const double dt)
{
    _lastContactCommand = _admittance.compute(_contactState, _contactReference, dt);
    return _admittance.apply_to_twist(baseTwist, _contactState, _lastContactCommand);
}

} } // namespace RobotLibrary::Control
