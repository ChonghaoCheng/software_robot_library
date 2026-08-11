/**
 * @file    AdmittanceContactController.h
 * @brief   Baseline admittance-style normal-force contact controller.
 */

#ifndef ADMITTANCE_CONTACT_CONTROLLER_H
#define ADMITTANCE_CONTACT_CONTROLLER_H

#include <Control/ContactControllerBase.h>

namespace RobotLibrary { namespace Control {

/**
 * @brief Simple force-to-normal-velocity baseline for contact tasks.
 *
 * The controller computes:
 *   forceError = targetNormalForce - measuredNormalForce
 *   vRelative  = forceVelocityGain * forceError / forceResponseGain
 *   vRobotN    = surfaceNormalVelocity + clamp(vRelative)
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
