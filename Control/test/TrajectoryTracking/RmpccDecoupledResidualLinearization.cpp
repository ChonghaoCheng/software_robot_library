#include "RmpccCostGeometry.h"
#include "RmpccPrediction.h"
#include "RmpccResidualLinearization.h"

#include <Eigen/Core>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace {

using State = RobotLibrary::Control::RmpccStateVector;
using Residuals = RobotLibrary::Control::RmpccDecoupledResiduals;

Eigen::Vector<double,6> curved_tangent(const double progress)
{
    const double s = std::clamp(progress, 0.0, 1.0);
    Eigen::Vector<double,6> tangent;
    tangent << 0.18 * std::cos(1.7 * s),
               0.12 * std::sin(2.1 * s) + 0.03,
               0.05 + 0.08 * s,
               0.21 + 0.06 * s,
              -0.17 + 0.04 * s,
               0.13 - 0.03 * s;
    return tangent;
}

Eigen::Matrix<double,12,1> stacked(const Residuals &residual)
{
    Eigen::Matrix<double,12,1> value;
    value << residual.contour, residual.lag;
    return value;
}

Eigen::Matrix<double,12,7> independent_jacobian(
    const State &state, const double step)
{
    Eigen::Matrix<double,12,7> jacobian;
    for(int column = 0; column < 7; ++column)
    {
        State plus = state;
        State minus = state;
        plus(column) += step;
        minus(column) -= step;
        double denominator = 2.0 * step;
        if(column == 6)
        {
            plus(6) = std::clamp(plus(6), 0.0, 1.0);
            minus(6) = std::clamp(minus(6), 0.0, 1.0);
            denominator = plus(6) - minus(6);
        }
        jacobian.col(column) =
            (stacked(RobotLibrary::Control::rmpcc_decoupled_residuals(
                 plus, 1e-10, curved_tangent))
             - stacked(RobotLibrary::Control::rmpcc_decoupled_residuals(
                 minus, 1e-10, curved_tangent))) / denominator;
    }
    return jacobian;
}

double objective(
    const State &state,
    const Eigen::Matrix<double,6,6> &contourWeight,
    const Eigen::Matrix<double,6,6> &lagWeight)
{
    const Residuals residual =
        RobotLibrary::Control::rmpcc_decoupled_residuals(
            state, 1e-10, curved_tangent);
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
    state << 0.055, -0.031, 0.027, 0.62, -0.41, 0.33, 0.43;
    State direction;
    direction << -0.18, 0.25, 0.14, 0.30, -0.27, 0.21, 0.39;
    direction.normalize();

    Eigen::Matrix<double,6,6> contourWeight =
        Eigen::Matrix<double,6,6>::Identity();
    contourWeight.diagonal() << 900.0, 1050.0, 820.0, 38.0, 46.0, 41.0;
    Eigen::Matrix<double,6,6> lagWeight =
        Eigen::Matrix<double,6,6>::Identity();
    lagWeight.diagonal() << 9.0, 12.0, 8.0, 5.0, 7.0, 6.0;

    constexpr double modelStep = 1e-6;
    constexpr double checkStep = 2e-7;
    const auto corrected =
        RobotLibrary::Control::rmpcc_linearize_decoupled_residuals(
            state, 1e-10, modelStep, curved_tangent);
    Eigen::Matrix<double,12,7> correctedJacobian;
    correctedJacobian.topRows<6>() = corrected.contourJacobian;
    correctedJacobian.bottomRows<6>() = corrected.lagJacobian;
    const Eigen::Matrix<double,12,7> referenceJacobian =
        independent_jacobian(state, checkStep);
    const double correctedResidualError =
        (correctedJacobian - referenceJacobian).norm()
        / std::max(1.0, referenceJacobian.norm());

    const Eigen::Matrix<double,6,6> errorMap =
        RobotLibrary::Control::rmpcc_decoupled_error_jacobian(
            state.head<6>(), modelStep);
    const Eigen::Vector3d positionTangent = curved_tangent(state(6)).head<3>();
    const Eigen::Matrix3d lag = RobotLibrary::Control::rmpcc_metric_projection<3>(
        positionTangent, Eigen::Matrix3d::Identity(), 1e-10);
    Eigen::Matrix<double,6,6> lagProjection = Eigen::Matrix<double,6,6>::Zero();
    lagProjection.block<3,3>(0,0) = lag;
    const Eigen::Matrix<double,6,6> contourProjection =
        Eigen::Matrix<double,6,6>::Identity() - lagProjection;
    Eigen::Matrix<double,12,7> frozenJacobian = Eigen::Matrix<double,12,7>::Zero();
    frozenJacobian.topLeftCorner<6,6>() = contourProjection * errorMap;
    frozenJacobian.bottomLeftCorner<6,6>() = lagProjection * errorMap;
    const double frozenResidualError =
        (frozenJacobian - referenceJacobian).norm()
        / std::max(1.0, referenceJacobian.norm());

    const State correctedGradient =
        2.0 * corrected.contourJacobian.transpose()
              * contourWeight * corrected.residual.contour
        + 2.0 * corrected.lagJacobian.transpose()
              * lagWeight * corrected.residual.lag;
    const State frozenGradient =
        2.0 * frozenJacobian.topRows<6>().transpose()
              * contourWeight * corrected.residual.contour
        + 2.0 * frozenJacobian.bottomRows<6>().transpose()
              * lagWeight * corrected.residual.lag;
    const double finiteDifference =
        (objective(state + checkStep * direction, contourWeight, lagWeight)
         - objective(state - checkStep * direction, contourWeight, lagWeight))
        / (2.0 * checkStep);
    const double correctedObjectiveError = relative_error(
        correctedGradient.dot(direction), finiteDifference);
    const double frozenObjectiveError = relative_error(
        frozenGradient.dot(direction), finiteDifference);

    std::cout << std::setprecision(12)
              << "decoupled frozen residual Jacobian relative error: "
              << frozenResidualError << '\n'
              << "decoupled corrected residual Jacobian relative error: "
              << correctedResidualError << '\n'
              << "decoupled frozen objective directional relative error: "
              << frozenObjectiveError << '\n'
              << "decoupled corrected objective directional relative error: "
              << correctedObjectiveError << '\n';
    if(correctedResidualError > 1e-5 or correctedObjectiveError > 1e-5)
    {
        return 1;
    }
    if(frozenResidualError <= 1e-5 or frozenObjectiveError <= 1e-5)
    {
        std::cerr << "Curved-path fixture does not expose the frozen projector error.\n";
        return 2;
    }

    for(const double progress :
        std::array<double,5>{0.0, 0.5 * modelStep, 0.4,
                             1.0 - 0.5 * modelStep, 1.0})
    {
        State boundaryState = state;
        boundaryState(6) = progress;
        const auto boundary =
            RobotLibrary::Control::rmpcc_linearize_decoupled_residuals(
                boundaryState, 1e-10, modelStep, curved_tangent);
        Eigen::Matrix<double,12,7> boundaryJacobian;
        boundaryJacobian.topRows<6>() = boundary.contourJacobian;
        boundaryJacobian.bottomRows<6>() = boundary.lagJacobian;
        const double residualError =
            (boundaryJacobian - independent_jacobian(boundaryState, checkStep)).norm()
            / std::max(1.0, independent_jacobian(boundaryState, checkStep).norm());

        State feasibleDirection = direction;
        if(progress == 0.0)
        {
            feasibleDirection(6) = std::abs(feasibleDirection(6));
        }
        else if(progress == 1.0)
        {
            feasibleDirection(6) = -std::abs(feasibleDirection(6));
        }
        const State gradient =
            2.0 * boundary.contourJacobian.transpose()
                  * contourWeight * boundary.residual.contour
            + 2.0 * boundary.lagJacobian.transpose()
                  * lagWeight * boundary.residual.lag;
        double objectiveDifference = 0.0;
        if(progress == 0.0 or progress == 1.0)
        {
            objectiveDifference =
                (objective(boundaryState + checkStep * feasibleDirection,
                           contourWeight, lagWeight)
                 - objective(boundaryState, contourWeight, lagWeight)) / checkStep;
        }
        else
        {
            objectiveDifference =
                (objective(boundaryState + checkStep * feasibleDirection,
                           contourWeight, lagWeight)
                 - objective(boundaryState - checkStep * feasibleDirection,
                             contourWeight, lagWeight)) / (2.0 * checkStep);
        }
        const double objectiveError = relative_error(
            gradient.dot(feasibleDirection), objectiveDifference);
        std::cout << "decoupled boundary s=" << progress
                  << " residual error=" << residualError
                  << " objective error=" << objectiveError << '\n';
        if(residualError > 1e-5 or objectiveError > 1e-5)
        {
            return 3;
        }
    }
    return 0;
}
