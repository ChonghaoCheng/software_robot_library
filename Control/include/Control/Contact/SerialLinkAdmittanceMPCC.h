/**
 * @file SerialLinkAdmittanceMPCC.h
 * @brief MPCC motion command followed by an admittance normal-force merge.
 */

#ifndef SERIAL_LINK_ADMITTANCE_MPCC_H
#define SERIAL_LINK_ADMITTANCE_MPCC_H

#include <Control/Contact/AdmittanceContactController.h>
#include <Control/TrajectoryTracking/SerialLinkMPCC.h>

namespace RobotLibrary { namespace Control {

/**
 * Comparator baseline: solve the unchanged MPCC, replace only its base-frame
 * linear normal component, then run the shared resolved-rate joint controller.
 */
class SerialLinkAdmittanceMPCC : public SerialLinkMPCC
{
    public:
        using SerialLinkMPCC::SerialLinkMPCC;

        void set_contact_state(const ContactState &state) { _contactState = state; }
        void set_contact_reference(const ContactReference &reference) { _contactReference = reference; }

        const ContactDiagnostics &contact_diagnostics() const
        {
            return _admittance.diagnostics();
        }

        const ContactCommand &last_contact_command() const { return _lastContactCommand; }

    protected:
        Eigen::Vector<double,6>
        postprocess_base_twist(const Eigen::Vector<double,6> &baseTwist,
                               double dt) override;

    private:
        ContactState _contactState;
        ContactReference _contactReference;
        ContactCommand _lastContactCommand;
        AdmittanceContactController _admittance;
};

} } // namespace RobotLibrary::Control

#endif
