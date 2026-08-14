#include <Control/RmpccCostGeometry.h>
#include <Control/RmpccPrediction.h>
#include <Control/RmpccResidualLinearization.h>

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <iostream>

namespace {

using State = RobotLibrary::Control::RmpccStateVector;

Eigen::Vector<double,6> tangent(const double progress)
{
    Eigen::Vector<double,6> value;
    value << 0.16 + 0.09 * progress,
             -0.07 + 0.04 * progress,
              0.11 - 0.03 * progress,
              0.31 + 0.13 * progress,
             -0.23 + 0.10 * progress,
              0.17 - 0.07 * progress;
    return value;
}

double objective(const State &state,
                 const Eigen::Matrix<double,6,6> &metric,
                 const Eigen::Matrix<double,6,6> &contourWeight,
                 const Eigen::Matrix<double,6,6> &lagWeight)
{
    const auto residual = RobotLibrary::Control::rmpcc_full_screw_residuals(
        state, metric, 1e-10, tangent);
    return (residual.contour.transpose() * contourWeight * residual.contour)(0)
           + (residual.lag.transpose() * lagWeight * residual.lag)(0);
}

double relative_error(const double value, const double reference)
{
    return std::abs(value - reference) / std::max(1e-12, std::abs(reference));
}

} // namespace

int main()
{
    State state;
    state << 0.05, -0.04, 0.03, 0.58, -0.43, 0.36, 0.41;
    State direction;
    direction << -0.17, 0.23, 0.11, 0.31, -0.29, 0.19, 0.37;
    direction.normalize();

    Eigen::Matrix<double,6,6> metric = Eigen::Matrix<double,6,6>::Identity();
    metric.diagonal() << 1.0, 1.2, 0.9, 0.04, 0.06, 0.05;
    Eigen::Matrix<double,6,6> contourWeight =
        Eigen::Matrix<double,6,6>::Identity();
    contourWeight.diagonal() << 800.0, 1100.0, 950.0, 35.0, 48.0, 41.0;
    Eigen::Matrix<double,6,6> lagWeight =
        Eigen::Matrix<double,6,6>::Identity();
    lagWeight.diagonal() << 8.0, 12.0, 9.0, 6.0, 11.0, 7.0;

    constexpr double residualStep = 1e-6;
    const auto full =
        RobotLibrary::Control::rmpcc_linearize_full_screw_residuals(
            state, metric, 1e-10, residualStep, tangent);
    const State fullGradient =
        2.0 * full.contourJacobian.transpose()
              * contourWeight * full.residual.contour
        + 2.0 * full.lagJacobian.transpose()
              * lagWeight * full.residual.lag;

    const Eigen::Vector<double,6> errorTangent =
        RobotLibrary::Control::rmpcc_error_coordinate_path_tangent(
            state.head<6>(), tangent(state(6)));
    const auto frozenProjection = RobotLibrary::Control::rmpcc_error_projection(
        errorTangent, metric,
        RobotLibrary::Control::RmpccLagGeometry::FullScrew, 1e-10);
    Eigen::Matrix<double,6,7> frozenContour = Eigen::Matrix<double,6,7>::Zero();
    Eigen::Matrix<double,6,7> frozenLag = Eigen::Matrix<double,6,7>::Zero();
    frozenContour.leftCols<6>() = frozenProjection.contour;
    frozenLag.leftCols<6>() = frozenProjection.lag;
    const State frozenGradient =
        2.0 * frozenContour.transpose()
              * contourWeight * full.residual.contour
        + 2.0 * frozenLag.transpose()
              * lagWeight * full.residual.lag;

    constexpr double objectiveStep = 2e-6;
    const double finiteDifference =
        (objective(state + objectiveStep * direction,
                   metric, contourWeight, lagWeight)
         - objective(state - objectiveStep * direction,
                     metric, contourWeight, lagWeight))
        / (2.0 * objectiveStep);
    const double fullDerivative = fullGradient.dot(direction);
    const double frozenDerivative = frozenGradient.dot(direction);
    const double fullError = relative_error(fullDerivative, finiteDifference);
    const double frozenError = relative_error(frozenDerivative, finiteDifference);

    std::cout << "objective directional derivative finite difference: "
              << finiteDifference << '\n'
              << "frozen-projector directional relative error: "
              << frozenError << '\n'
              << "full-residual directional relative error: "
              << fullError << '\n';
    if(fullError > 1e-5)
    {
        std::cerr << "Full residual objective gradient failed directional check.\n";
        return 1;
    }
    if(frozenError <= 1e-5)
    {
        std::cerr << "Fixture does not distinguish frozen and full gradients.\n";
        return 2;
    }
    return 0;
}
