#include "RmpccCostGeometry.h"
#include "RmpccPrediction.h"
#include "RmpccResidualLinearization.h"

#include <Eigen/Core>

#include <array>
#include <algorithm>
#include <iomanip>
#include <iostream>

namespace {

using State = Eigen::Matrix<double,7,1>;
using Residual = Eigen::Matrix<double,12,1>;

Eigen::Matrix<double,6,1> tangent(const double progress)
{
    const double s = std::clamp(progress, 0.0, 1.0);
    Eigen::Matrix<double,6,1> value;
    value << 0.18 + 0.07 * s,
             -0.05 + 0.03 * s,
              0.09 - 0.02 * s,
              0.24 + 0.11 * s,
             -0.19 + 0.08 * s,
              0.13 - 0.05 * s;
    return value;
}

Residual residual(
    const State &state,
    const Eigen::Matrix<double,6,6> &metric,
    const RobotLibrary::Control::RmpccLagGeometry geometry =
        RobotLibrary::Control::RmpccLagGeometry::FullScrew)
{
    const Eigen::Matrix<double,6,1> g =
        RobotLibrary::Control::rmpcc_error_coordinate_path_tangent(
            state.head<6>(), tangent(state(6)));
    const auto projection = RobotLibrary::Control::rmpcc_error_projection(
        g, metric, geometry, 1e-10);
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

    for(const double progress : std::array<double,4>{0.0, 1.0, 1.0 - 0.5 * h, 0.4})
    {
        State boundaryState = state;
        boundaryState(6) = progress;
        Eigen::Matrix<double,12,7> feasibleDifference;
        for(int column = 0; column < 7; ++column)
        {
            State plus = boundaryState;
            State minus = boundaryState;
            plus(column) += h;
            minus(column) -= h;
            double denominator = 2.0 * h;
            if(column == 6)
            {
                plus(6) = std::clamp(plus(6), 0.0, 1.0);
                minus(6) = std::clamp(minus(6), 0.0, 1.0);
                denominator = plus(6) - minus(6);
            }
            feasibleDifference.col(column) =
                (residual(plus, metric) - residual(minus, metric)) / denominator;
        }

        const auto boundaryLinearization =
            RobotLibrary::Control::rmpcc_linearize_full_screw_residuals(
                boundaryState, metric, 1e-10, h, tangent);
        Eigen::Matrix<double,12,7> boundaryJacobian;
        boundaryJacobian.topRows<6>() = boundaryLinearization.contourJacobian;
        boundaryJacobian.bottomRows<6>() = boundaryLinearization.lagJacobian;
        const double boundaryError =
            (boundaryJacobian - feasibleDifference).norm()
            / std::max(1.0, feasibleDifference.norm());
        std::cout << "boundary residual Jacobian s=" << std::setprecision(10) << progress
                  << " relative error: " << boundaryError << '\n';
        if(boundaryError > 1e-5)
        {
            std::cerr << "Full residual boundary Jacobian mismatch at s="
                      << progress << ": relative error=" << boundaryError << '\n';
            return 3;
        }
    }

    using RobotLibrary::Control::RmpccLagGeometry;
    for(const RmpccLagGeometry geometry :
        std::array<RmpccLagGeometry,3>{
            RmpccLagGeometry::SplitTranslationRotation,
            RmpccLagGeometry::TranslationOnly,
            RmpccLagGeometry::RotationOnly})
    {
        for(const double progress :
            std::array<double,4>{0.0, 1.0, 1.0 - 0.5 * h, 0.4})
        {
            State geometryState = state;
            geometryState(6) = progress;
            Eigen::Matrix<double,12,7> expected;
            for(int column = 0; column < 7; ++column)
            {
                State plus = geometryState;
                State minus = geometryState;
                plus(column) += h;
                minus(column) -= h;
                double denominator = 2.0 * h;
                if(column == 6)
                {
                    plus(6) = std::clamp(plus(6), 0.0, 1.0);
                    minus(6) = std::clamp(minus(6), 0.0, 1.0);
                    denominator = plus(6) - minus(6);
                }
                expected.col(column) =
                    (residual(plus, metric, geometry)
                     - residual(minus, metric, geometry)) / denominator;
            }
            const auto linearization =
                RobotLibrary::Control::rmpcc_linearize_projected_residuals(
                    geometryState, metric, geometry, 1e-10, h, tangent);
            Eigen::Matrix<double,12,7> actual;
            actual.topRows<6>() = linearization.contourJacobian;
            actual.bottomRows<6>() = linearization.lagJacobian;
            const double error = (actual - expected).norm()
                / std::max(1.0, expected.norm());
            if(error > 1e-5)
            {
                std::cerr << "Projected residual Jacobian mismatch for geometry "
                          << static_cast<int>(geometry) << " at s=" << progress
                          << ": " << error << '\n';
                return 4;
            }
        }
    }
    return 0;
}
