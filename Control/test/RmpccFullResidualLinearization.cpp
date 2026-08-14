#include <Control/RmpccCostGeometry.h>
#include <Control/RmpccPrediction.h>
#include <Control/RmpccResidualLinearization.h>

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
    Eigen::Matrix<double,12,7> frozenJacobian = Eigen::Matrix<double,12,7>::Zero();
    frozenJacobian.block<6,6>(0,0) = frozen.contour;
    frozenJacobian.block<6,6>(6,0) = frozen.lag;

    const auto full =
        RobotLibrary::Control::rmpcc_linearize_full_screw_residuals(
            state, metric, 1e-10, h, tangent);
    Eigen::Matrix<double,12,7> fullJacobian;
    fullJacobian.topRows<6>() = full.contourJacobian;
    fullJacobian.bottomRows<6>() = full.lagJacobian;

    const double frozenRelativeError =
        (frozenJacobian - finiteDifference).norm()
        / std::max(1.0, finiteDifference.norm());
    const double fullRelativeError =
        (fullJacobian - finiteDifference).norm()
        / std::max(1.0, finiteDifference.norm());
    std::cout << "frozen-projector residual gradient relative error: "
              << frozenRelativeError << '\n'
              << "full-residual gradient relative error: "
              << fullRelativeError << '\n';
    if(fullRelativeError > 1e-5)
    {
        std::cerr << "Complete residual Jacobian does not match finite difference.\n";
        return 1;
    }
    if(frozenRelativeError <= 1e-5)
    {
        std::cerr << "Fixture no longer exposes the frozen-projector approximation.\n";
        return 2;
    }
    return 0;
}
