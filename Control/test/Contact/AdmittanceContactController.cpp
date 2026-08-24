/**
 * @file AdmittanceContactController.cpp
 * @brief Deterministic sign, feedforward, and twist-merge tests.
 */

#include <Control/Contact/AdmittanceContactController.h>

#include <cmath>
#include <iostream>
#include <string>

namespace {

int failures = 0;

void check(const bool condition, const std::string &message)
{
    if(!condition)
    {
        std::cerr << "FAIL: " << message << "\n";
        ++failures;
    }
}

bool near(const double lhs, const double rhs, const double tolerance = 1e-12)
{
    return std::abs(lhs - rhs) <= tolerance;
}

} // namespace

int main()
{
    using namespace RobotLibrary::Control;

    AdmittanceContactController controller;
    ContactState state;
    state.normalBase = Eigen::Vector3d::UnitZ();
    state.inContact = true;

    ContactReference reference;
    reference.targetNormalForce = 5.0;
    reference.forceResponseGain = 1000.0;
    reference.forceVelocityGain = 1.0;
    reference.maxNormalVelocity = 1.0;

    state.measuredNormalForce = 3.0;
    ContactCommand command = controller.compute(state, reference, 0.001);
    check(command.desiredRelativeNormalVelocity > 0.0,
          "low compressive force commands motion along the inward normal");

    state.measuredNormalForce = 7.0;
    command = controller.compute(state, reference, 0.001);
    check(command.desiredRelativeNormalVelocity < 0.0,
          "high compressive force commands motion opposite the inward normal");

    state.measuredNormalForce = reference.targetNormalForce;
    state.surfaceVelocityBase = 0.01 * Eigen::Vector3d::UnitZ();
    command = controller.compute(state, reference, 0.001);
    check(near(command.desiredRelativeNormalVelocity, 0.0),
          "zero force error produces zero board-relative velocity");
    check(near(command.desiredRobotNormalVelocity, 0.01),
          "surface normal velocity is passed through as robot feedforward");

    Eigen::Vector<double,6> motionTwist;
    motionTwist << 0.12, -0.07, 0.31, 0.4, -0.5, 0.6;
    command.active = true;
    command.desiredRobotNormalVelocity = -0.02;
    const Eigen::Vector<double,6> merged =
        controller.apply_to_twist(motionTwist, state, command);

    check(near(merged.x(), motionTwist.x()) && near(merged.y(), motionTwist.y()),
          "twist merge preserves both tangent linear components");
    check((merged.tail<3>() - motionTwist.tail<3>()).norm() <= 1e-12,
          "twist merge preserves every angular component");
    check(near(merged.z(), command.desiredRobotNormalVelocity),
          "twist merge replaces only the linear normal component");

    if(failures == 0)
    {
        std::cout << "admittance_contact_controller_test PASS\n";
        return 0;
    }
    std::cerr << "admittance_contact_controller_test FAIL (" << failures << " checks)\n";
    return 1;
}
