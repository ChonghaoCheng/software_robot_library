#include <Control/RmpccCostGeometry.h>
#include <Control/RmpccPrediction.h>

#include <Eigen/Core>
#include <iostream>

namespace {

using State = Eigen::Matrix<double,7,1>;
using Residual = Eigen::Matrix<double,12,1>;

Eigen::Matrix<double,6,1> tangent(const double progress)
{
    Eigen::Matrix<double,6,1> value;
    value << 0.18 + 0.07 * progress,
             -0.05 + 0.03 * progress,
              0.09 - 0.02 * progress,
              0.24 + 0.11 * progress,
             -0.19 + 0.08 * progress,
              0.13 - 0.05 * progress;
    return value;
}

Residual residual(const State &state, const Eigen::Matrix<double,6,6> &metric)
{
    const Eigen::Matrix<double,6,1> g =
        RobotLibrary::Control::rmpcc_error_coordinate_path_tangent(
            state.head<6>(), tangent(state(6)));
    const auto projection = RobotLibrary::Control::rmpcc_error_projection(
        g, metric, RobotLibrary::Control::RmpccLagGeometry::FullScrew, 1e-10);
    Residual value;
    value.head<6>() = projection.contour * state.head<6>();
    value.tail<6>() = projection.lag * state.head<6>();
    return value;
}

} // namespace

int main()
{
    State state;
    state << 0.04, -0.03, 0.02, 0.42, -0.31, 0.27, 0.46;
    Eigen::Matrix<double,6,6> metric = Eigen::Matrix<double,6,6>::Identity();
    metric.diagonal() << 1.0, 1.3, 0.8, 0.04, 0.06, 0.05;

    constexpr double h = 1e-6;
    Eigen::Matrix<double,12,7> finiteDifference;
    for(int column = 0; column < 7; ++column)
    {
        State plus = state;
        State minus = state;
        plus(column) += h;
        minus(column) -= h;
        finiteDifference.col(column) =
            (residual(plus, metric) - residual(minus, metric)) / (2.0 * h);
    }

    const Eigen::Matrix<double,6,1> g =
        RobotLibrary::Control::rmpcc_error_coordinate_path_tangent(
            state.head<6>(), tangent(state(6)));
    const auto frozen = RobotLibrary::Control::rmpcc_error_projection(
        g, metric, RobotLibrary::Control::RmpccLagGeometry::FullScrew, 1e-10);
    Eigen::Matrix<double,12,7> implemented = Eigen::Matrix<double,12,7>::Zero();
    implemented.block<6,6>(0,0) = frozen.contour;
    implemented.block<6,6>(6,0) = frozen.lag;

    const double relativeError =
        (implemented - finiteDifference).norm()
        / std::max(1.0, finiteDifference.norm());
    std::cout << "full-screw residual gradient relative error: "
              << relativeError << '\n';
    if(relativeError > 1e-5)
    {
        std::cerr << "Frozen projector omits state-dependent residual derivatives.\n";
        return 1;
    }
    return 0;
}
