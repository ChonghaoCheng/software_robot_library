#include <Math/DiscreteIntegratorLQR.h>

#include <cmath>

int main()
{
    constexpr double dt = 0.002;
    constexpr double gain = 20.0;
    constexpr double controlWeight = 0.1;
    const double q = RobotLibrary::Math::integrator_stage_weight_for_gain(
        gain, controlWeight, dt);
    const double terminal = RobotLibrary::Math::integrator_terminal_weight_for_gain(
        gain, controlWeight, dt);
    const double riccati = RobotLibrary::Math::integrator_terminal_weight(
        q, controlWeight, dt);
    if(std::abs(terminal - riccati) > 1e-10) return 1;

    const double recoveredGain =
        dt * terminal / (controlWeight + dt * dt * terminal);
    if(std::abs(recoveredGain - gain) > 1e-10) return 2;
    return 0;
}
