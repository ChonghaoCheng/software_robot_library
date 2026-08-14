#include <Eigen/Core>
#include <iostream>

int main()
{
    constexpr int horizon = 3;
    constexpr double dt = 0.1;
    constexpr double remaining = 1.0;
    constexpr double scheduleRemaining = 0.05;
    const Eigen::Vector3d progressRates(0.5, 0.5, 0.5);

    const double completionAdvance = dt * progressRates.sum();
    const double firstStepAdvance = dt * progressRates(0);
    const double implementedScheduleAdvance = dt * progressRates.sum();

    if(completionAdvance > remaining + 1e-12)
    {
        std::cerr << "Test fixture violates the full-horizon completion constraint.\n";
        return 2;
    }
    if(firstStepAdvance > scheduleRemaining + 1e-12)
    {
        std::cerr << "Test fixture violates documented first-step schedule semantics.\n";
        return 3;
    }
    if(implementedScheduleAdvance > scheduleRemaining + 1e-12)
    {
        std::cerr << "Current schedule row constrains the full horizon ("
                  << implementedScheduleAdvance << ") instead of only sdot_0 ("
                  << firstStepAdvance << ").\n";
        return 1;
    }
    return 0;
}
