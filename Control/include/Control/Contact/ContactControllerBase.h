/**
 * @file    ContactControllerBase.h
 * @brief   Base class for reusable normal-force contact-control baselines.
 */

#ifndef CONTACT_CONTROLLER_BASE_H
#define CONTACT_CONTROLLER_BASE_H

#include <Control/Contact/ContactDataStructures.h>

#include <Eigen/Core>

namespace RobotLibrary { namespace Control {

/**
 * @brief Base interface for contact controllers that can be composed with motion controllers.
 *
 * This class intentionally has no ROS, TF, logger, or sensor subscriptions. A
 * front-end supplies ContactState and ContactReference; the controller returns a
 * normal-velocity command that can be merged with any motion-controller twist.
 */
class ContactControllerBase
{
    public:
        virtual ~ContactControllerBase() = default;

        /**
         * @brief Compute the contact command for the current contact state.
         */
        virtual ContactCommand
        compute(const ContactState &state,
                const ContactReference &reference,
                double dt) = 0;

        /**
         * @brief Merge a contact command into a 6D motion twist by replacing only
         *        the linear normal component and preserving tangent/angular motion.
         */
        Eigen::Vector<double,6>
        apply_to_twist(const Eigen::Vector<double,6> &motionTwist,
                       const ContactState &state,
                       const ContactCommand &command);

        /**
         * @brief Reset controller memory and diagnostics.
         */
        virtual void
        reset();

        /**
         * @brief Latest diagnostics from compute() or apply_to_twist().
         */
        const ContactDiagnostics&
        diagnostics() const { return _diagnostics; }

    protected:
        ContactDiagnostics _diagnostics;

        static bool
        normal_is_valid(const Eigen::Vector3d &normal);

        static Eigen::Vector3d
        unit_normal(const Eigen::Vector3d &normal);
};

} } // namespace

#endif
