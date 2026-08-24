/** @file ContactClosedLoop.cpp
 * Deterministic compliant-contact plant tests with mismatch, delay, noise, and re-contact.
 */

#include <Control/Contact/AdmittanceContactController.h>

#include <algorithm>
#include <cmath>
#include <deque>
#include <iostream>
#include <string>

namespace {

using namespace RobotLibrary::Control;

int failures = 0;

void check(const bool condition, const std::string &message)
{
    if(!condition)
    {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

double compressive_force(const double stiffness,
                         const double toolCoordinate,
                         const double surfaceCoordinate)
{
    return stiffness * std::max(0.0, toolCoordinate - surfaceCoordinate);
}

} // namespace

int main()
{
    constexpr double pi = 3.14159265358979323846;
    constexpr double dt = 0.001;

    // Plant stiffness differs from the controller model, force feedback is
    // delayed by 10 ms, and deterministic sensor noise is present.
    {
        AdmittanceContactController controller;
        ContactReference reference;
        reference.targetNormalForce = 5.0;
        reference.forceResponseGain = 1500.0;
        reference.forceVelocityGain = 4.0;
        reference.maxNormalVelocity = 0.03;
        reference.forceDeadband = 0.02;

        constexpr double trueStiffness = 800.0;
        double surface = 0.0;
        double tool = 0.002;
        double trueForce = compressive_force(trueStiffness, tool, surface);
        std::deque<double> delayedForce(10, trueForce);
        double finalWindowAbsoluteError = 0.0;
        int finalWindowSamples = 0;
        bool velocityBoundRespected = true;

        for(int step = 0; step < 4000; ++step)
        {
            const double time = step * dt;
            const double surfaceVelocity = 0.01 * std::sin(2.0 * pi * 0.5 * time);
            trueForce = compressive_force(trueStiffness, tool, surface);
            delayedForce.push_back(trueForce);
            const double delayedMeasurement = delayedForce.front();
            delayedForce.pop_front();
            const double measuredForce = std::max(
                0.0, delayedMeasurement + 0.05 * std::sin(2.0 * pi * 17.0 * time));

            ContactState state;
            state.normalBase = Eigen::Vector3d::UnitZ();
            state.surfaceVelocityBase = surfaceVelocity * Eigen::Vector3d::UnitZ();
            state.measuredNormalForce = measuredForce;
            state.inContact = true;
            const ContactCommand command = controller.compute(state, reference, dt);
            velocityBoundRespected = velocityBoundRespected
                && std::abs(command.desiredRelativeNormalVelocity)
                    <= reference.maxNormalVelocity + 1e-12;

            tool += dt * command.desiredRobotNormalVelocity;
            surface += dt * surfaceVelocity;
            if(step >= 3500)
            {
                finalWindowAbsoluteError +=
                    std::abs(reference.targetNormalForce - trueForce);
                ++finalWindowSamples;
            }
        }

        const double meanFinalError =
            finalWindowAbsoluteError / static_cast<double>(finalWindowSamples);
        check(velocityBoundRespected,
              "normal relative velocity remains bounded under delay and noise");
        check(meanFinalError < 0.12,
              "force converges with plant/model stiffness mismatch, delay, and noise");
    }

    // Start outside contact, approach, establish contact, then regulate force.
    {
        AdmittanceContactController controller;
        ContactReference reference;
        reference.targetNormalForce = 3.0;
        reference.forceResponseGain = 1500.0;
        reference.forceVelocityGain = 4.0;
        reference.maxNormalVelocity = 0.03;
        reference.approachWhenNotInContact = true;
        reference.approachVelocity = 0.01;

        constexpr double trueStiffness = 1000.0;
        double surface = 0.0;
        double tool = -0.002;
        bool contactEstablished = false;
        bool approachWasActive = false;
        double force = 0.0;

        for(int step = 0; step < 2500; ++step)
        {
            force = compressive_force(trueStiffness, tool, surface);
            ContactState state;
            state.normalBase = Eigen::Vector3d::UnitZ();
            state.measuredNormalForce = force;
            state.inContact = force > 1e-9;
            contactEstablished = contactEstablished || state.inContact;
            const ContactCommand command = controller.compute(state, reference, dt);
            approachWasActive = approachWasActive || command.approachActive;
            tool += dt * command.desiredRobotNormalVelocity;
        }

        check(approachWasActive, "free-space approach phase is exercised");
        check(contactEstablished, "approach reaches the compliant surface");
        check(std::abs(force - reference.targetNormalForce) < 0.08,
              "controller regulates force after re-contact");
    }

    if(failures == 0)
    {
        std::cout << "contact_closed_loop_test PASS\n";
        return 0;
    }
    std::cerr << "contact_closed_loop_test FAIL (" << failures << " checks)\n";
    return 1;
}
