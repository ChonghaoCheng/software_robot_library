#include "ProgressSchedule.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>

int main()
{
    using RobotLibrary::Control::reference_schedule_progress;
    using RobotLibrary::Control::reference_schedule_rates;
    constexpr double dt = 0.002;
    constexpr double duration = 1.003;

    if(reference_schedule_progress(0, dt, duration) != 0.0
       or reference_schedule_progress(10000, dt, duration) != 1.0)
    {
        std::cerr << "Reference schedule endpoint semantics are wrong.\n";
        return 1;
    }

    for(const std::uint64_t step : {0ULL, 17ULL, 500ULL, 501ULL, 502ULL})
    {
        const double progress = reference_schedule_progress(step, dt, duration);
        const Eigen::VectorXd rates = reference_schedule_rates(
            progress, 7, dt, duration);
        double integrated = progress;
        for(int stage = 0; stage < rates.size(); ++stage)
        {
            integrated += dt * rates(stage);
            const double expected = std::clamp(
                progress + static_cast<double>(stage + 1) * dt / duration,
                0.0, 1.0);
            if(std::abs(integrated - expected) > 1e-14)
            {
                std::cerr << "Horizon rate does not reproduce supplied schedule.\n";
                return 2;
            }
        }
    }
    return 0;
}
