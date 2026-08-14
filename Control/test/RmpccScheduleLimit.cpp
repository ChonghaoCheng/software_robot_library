#include <Control/RmpccProgressConstraints.h>
#include <Eigen/Core>
#include <iostream>

int main()
{
    constexpr int horizon = 3;
    constexpr double dt = 0.1;
    constexpr double remaining = 1.0;
    constexpr double scheduleRemaining = 0.05;
    const Eigen::Vector3d progressRates(0.5, 0.5, 0.5);

    const Eigen::VectorXd completionRow =
        RobotLibrary::Control::rmpcc_completion_progress_row(horizon, dt);
    const Eigen::VectorXd scheduleRow =
        RobotLibrary::Control::rmpcc_schedule_progress_row(horizon, dt);
    const double completionAdvance = completionRow.dot(progressRates);
    const double firstStepAdvance = dt * progressRates(0);
    const double implementedScheduleAdvance = scheduleRow.dot(progressRates);

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
        std::cerr << "Schedule row does not implement first-step semantics.\n";
        return 1;
    }

    Eigen::VectorXd clipped = progressRates;
    const Eigen::VectorXd lower = Eigen::Vector3d::Zero();
    const Eigen::VectorXd upper = Eigen::Vector3d::Ones();
    RobotLibrary::Control::rmpcc_clip_progress_rates(
        clipped, lower, upper, dt, remaining, scheduleRemaining);
    if(std::abs(clipped(0) - 0.5) > 1e-12
       or std::abs(clipped.tail(2).sum() - 1.0) > 1e-12)
    {
        std::cerr << "Warm-start clipping still applies scheduleRemaining to the horizon.\n";
        return 4;
    }
    return 0;
}
