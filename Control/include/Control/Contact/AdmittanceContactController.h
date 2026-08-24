/**
 * @file    AdmittanceContactController.h
 * @brief   Baseline admittance-style normal-force contact controller.
 */

#ifndef ADMITTANCE_CONTACT_CONTROLLER_H
#define ADMITTANCE_CONTACT_CONTROLLER_H

#include <Control/Contact/ContactControllerBase.h>

namespace RobotLibrary { namespace Control {

/**
 * @brief Simple force-to-normal-velocity baseline for contact tasks.
 *
 * The controller computes:
 *   forceError = targetNormalForce - measuredNormalForce
 *   vRelative  = forceVelocityGain * forceError / forceResponseGain
 *   vRobotN    = surfaceNormalVelocity + clamp(vRelative)
 *
 * Positive normal velocity increases compressive force, and measured force is
 * positive in compression. apply_to_twist() replaces only the linear normal
 * component; tangent linear and all angular components are preserved.
 *
 * It is intended as a baseline module that MPC, MPCC, RMPCC, PID, or future
 * hybrid motion-force controllers can call by composition.
 */
class AdmittanceContactController : public ContactControllerBase
{
    public:
        ContactCommand
        compute(const ContactState &state,
                const ContactReference &reference,
                double dt) override;

    private:
        static double
        apply_deadband(double value, double deadband);
};

} } // namespace

#endif
