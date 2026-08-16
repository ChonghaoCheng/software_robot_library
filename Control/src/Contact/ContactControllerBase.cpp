/**
 * @file    ContactControllerBase.cpp
 * @brief   Base utilities for reusable contact-control baselines.
 */

#include <Control/Contact/ContactControllerBase.h>

#include <cmath>

namespace RobotLibrary { namespace Control {

bool
ContactControllerBase::normal_is_valid(const Eigen::Vector3d &normal)
{
    return normal.allFinite() && normal.norm() > 1e-9;
}

Eigen::Vector3d
ContactControllerBase::unit_normal(const Eigen::Vector3d &normal)
{
    if(!normal_is_valid(normal))
    {
        return Eigen::Vector3d::Zero();
    }

    return normal.normalized();
}

Eigen::Vector<double,6>
ContactControllerBase::apply_to_twist(const Eigen::Vector<double,6> &motionTwist,
                                      const ContactState &state,
                                      const ContactCommand &command)
{
    if(!command.active || !normal_is_valid(state.normalBase))
    {
        return motionTwist;
    }

    const Eigen::Vector3d normal = unit_normal(state.normalBase);
    const Eigen::Vector3d linearMotion = motionTwist.head<3>();
    const double motionNormalBefore = normal.dot(linearMotion);
    const Eigen::Vector3d tangentMotion = linearMotion - normal * motionNormalBefore;

    Eigen::Vector<double,6> mergedTwist = motionTwist;
    mergedTwist.head<3>() = tangentMotion + normal * command.desiredRobotNormalVelocity;

    _diagnostics.motionNormalVelocityBeforeMerge = motionNormalBefore;
    _diagnostics.motionNormalVelocityAfterMerge = normal.dot(mergedTwist.head<3>());

    return mergedTwist;
}

void
ContactControllerBase::reset()
{
    _diagnostics = ContactDiagnostics{};
}

} } // namespace
