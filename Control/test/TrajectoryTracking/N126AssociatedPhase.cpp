#include "RmpccAssociatedPhase.h"
#include "RmpccPhaseResidual.h"

#include <Math/MathFunctions.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

namespace {

using namespace RobotLibrary::Control;
constexpr double kPi = 3.14159265358979323846;

Eigen::Matrix4d reference(const double s)
{
    Eigen::Matrix4d result = Eigen::Matrix4d::Identity();
    result.block<3,1>(0,3) << 0.25 * s,
        0.04 * std::sin(2.0 * kPi * s),
        0.03 * std::sin(kPi * s);
    const Eigen::Vector3d axis = Eigen::Vector3d(0.3, 0.8, 0.5).normalized();
    result.block<3,3>(0,0) =
        (Eigen::AngleAxisd(0.9 * std::sin(kPi * s), axis)
         * Eigen::AngleAxisd(0.5 * s, Eigen::Vector3d::UnitZ()))
            .toRotationMatrix();
    Eigen::Matrix4d attachment = Eigen::Matrix4d::Identity();
    attachment(2,3) = 0.05;
    return result * attachment;
}

Eigen::Vector<double,6> tangent(const double s)
{
    constexpr double h = 1e-6;
    const double lo = std::max(0.0, s - h);
    const double hi = std::min(1.0, s + h);
    return RobotLibrary::Math::se3_logarithm(
        RobotLibrary::Math::se3_inverse(reference(lo)) * reference(hi))
        / (hi - lo);
}

Eigen::Matrix<double,6,6> lag_weight()
{
    Eigen::Matrix<double,6,6> result =
        Eigen::Matrix<double,6,6>::Zero();
    result.diagonal() << 10.0, 10.0, 10.0, 0.4, 0.4, 0.4;
    return result;
}

double fitted_order(const std::vector<double> &scales,
                    const std::vector<double> &errors)
{
    double meanX = 0.0;
    double meanY = 0.0;
    for(std::size_t i = 0; i < scales.size(); ++i)
    {
        meanX += std::log(scales[i]);
        meanY += std::log(std::max(errors[i], 1e-18));
    }
    meanX /= scales.size();
    meanY /= scales.size();
    double numerator = 0.0;
    double denominator = 0.0;
    for(std::size_t i = 0; i < scales.size(); ++i)
    {
        const double x = std::log(scales[i]) - meanX;
        numerator += x * (std::log(std::max(errors[i], 1e-18)) - meanY);
        denominator += x * x;
    }
    return numerator / denominator;
}

void require(const bool condition, const char *message)
{
    if(not condition) throw std::runtime_error(message);
}

} // namespace

int main()
{
    try
    {
        const auto transform = [](const double s) { return reference(s); };
        const auto derivative = [](const double s) { return tangent(s); };
        RmpccPoseArcTable arc;
        RmpccPoseArcTable referenceArc;
        arc.build(lag_weight(), derivative, 4097);
        referenceArc.build(lag_weight(), derivative, 16385);

        std::mt19937_64 generator(12620260817ULL);
        std::uniform_real_distribution<double> progress(0.0, 1.0);
        double maximumArcDifference = 0.0;
        for(int sample = 0; sample < 200; ++sample)
        {
            const double a = progress(generator);
            const double b = progress(generator);
            maximumArcDifference = std::max(maximumArcDifference,
                std::abs(arc.lag(a, b) - referenceArc.lag(a, b)));
        }
        require(maximumArcDifference < 1e-7, "pose-arc table accuracy failed");

        const std::vector<double> scales{0.008, 0.004, 0.002, 0.001};
        std::vector<double> phaseErrors;
        std::vector<double> exactContours;
        std::vector<double> localExactErrors;
        std::vector<double> arcCostErrors;
        std::vector<double> unifiedDecoupledErrors;
        for(const double scale : scales)
        {
            double phaseSquared = 0.0;
            double contourSquared = 0.0;
            double localExactSquared = 0.0;
            double arcSquared = 0.0;
            double decoupledSquared = 0.0;
            int count = 0;
            for(int point = 0; point < 25; ++point)
            {
                const double s = 0.15 + 0.70 * point / 24.0;
                for(const double sign : {-1.0, 1.0})
                {
                    const double ds = sign * scale;
                    RmpccStateVector state = RmpccStateVector::Zero();
                    state(6) = s;
                    state.head<6>() = RobotLibrary::Math::se3_logarithm(
                        RobotLibrary::Math::se3_inverse(reference(s))
                        * reference(s + ds));
                    const auto exact = rmpcc_associated_residuals(
                        state,
                        RmpccContourResidualGeometry::ExactAssociatedUnifiedSE3,
                        1e-12, arc, transform, derivative);
                    const auto decoupled = rmpcc_associated_residuals(
                        state,
                        RmpccContourResidualGeometry::AssociatedDecoupledCartesianSO3,
                        1e-12, arc, transform, derivative);
                    const auto local = rmpcc_phase_residuals(
                        state, Eigen::Matrix<double,6,6>::Identity(),
                        RmpccPhaseAssociation::TaskPointXYZ, 1e-8, 1e-12,
                        transform, derivative);
                    phaseSquared += std::pow(
                        exact.context.phaseCorrection - ds, 2);
                    contourSquared += exact.contour.squaredNorm();
                    localExactSquared +=
                        (exact.contour - local.contour).squaredNorm();
                    const Eigen::Vector<double,6> g =
                        rmpcc_error_coordinate_path_tangent(
                            state.head<6>(), tangent(s));
                    const double localCost = ds * ds
                        * (g.transpose() * lag_weight() * g)(0);
                    const double arcCost = arc.q_lag()
                        * exact.scalarPoseArcLag * exact.scalarPoseArcLag;
                    arcSquared += std::pow(arcCost - localCost, 2);
                    decoupledSquared +=
                        (decoupled.contour - exact.contour).squaredNorm();

                    const Eigen::Vector<double,6> localLagVector =
                        g * exact.context.phaseCorrection;
                    const double localPoseLag = std::copysign(
                        std::sqrt((localLagVector.transpose() * lag_weight()
                                   * localLagVector)(0) / arc.q_lag()),
                        exact.context.phaseCorrection);
                    const double vectorCost =
                        ((g * exact.context.phaseCorrection).transpose()
                         * lag_weight()
                         * (g * exact.context.phaseCorrection))(0);
                    require(std::abs(arc.q_lag() * localPoseLag * localPoseLag
                                     - vectorCost) < 1e-12,
                            "local scalar/vector identity failed");
                    ++count;
                }
            }
            phaseErrors.push_back(std::sqrt(phaseSquared / count));
            exactContours.push_back(std::sqrt(contourSquared / count));
            localExactErrors.push_back(std::sqrt(localExactSquared / count));
            arcCostErrors.push_back(std::sqrt(arcSquared / count));
            unifiedDecoupledErrors.push_back(
                std::sqrt(decoupledSquared / count));
        }

        const double phaseOrder = fitted_order(scales, phaseErrors);
        const double contourOrder = fitted_order(scales, exactContours);
        const double localExactOrder = fitted_order(scales, localExactErrors);
        const double arcOrder = fitted_order(scales, arcCostErrors);
        const double decoupledOrder = fitted_order(scales, unifiedDecoupledErrors);
        require(phaseOrder >= 1.8, "pure-progress phase order failed");
        require(contourOrder >= 1.8, "pure-progress exact contour order failed");
        require(localExactOrder >= 1.8, "local/exact contour order failed");
        require(arcOrder >= 2.8, "pose-arc local equivalence order failed");
        require(decoupledOrder >= 1.8, "unified/decoupled local order failed");

        RmpccStateVector state = RmpccStateVector::Zero();
        state(6) = 0.43;
        state.head<6>() << 0.002, -0.001, 0.0015, 0.015, -0.01, 0.02;
        const auto tpa = rmpcc_associated_residuals(
            state, RmpccContourResidualGeometry::ExactAssociatedUnifiedSE3,
            1e-12, arc, transform, derivative);
        const auto dt = rmpcc_associated_residuals(
            state, RmpccContourResidualGeometry::AssociatedDecoupledCartesianSO3,
            1e-12, arc, transform, derivative);
        require(tpa.context.internalProgress == dt.context.internalProgress
                && tpa.context.phaseCorrection == dt.context.phaseCorrection
                && tpa.context.associatedProgress == dt.context.associatedProgress
                && (tpa.context.referenceAtAssociatedProgress
                    - dt.context.referenceAtAssociatedProgress).norm() == 0.0
                && tpa.scalarPoseArcLag == dt.scalarPoseArcLag
                && tpa.context.observable == dt.context.observable,
                "TPA/DT shared context identity failed");

        double maximumDirectionalDerivativeError = 0.0;
        for(const auto geometry : {
                RmpccContourResidualGeometry::ExactAssociatedUnifiedSE3,
                RmpccContourResidualGeometry::AssociatedDecoupledCartesianSO3})
        {
            const auto linearized = rmpcc_linearize_associated_residuals(
                state, geometry, 1e-12, 1e-6, arc, transform, derivative);
            Eigen::Vector<double,7> direction;
            direction << 0.2, -0.3, 0.1, 0.4, -0.2, 0.15, 0.25;
            direction.normalize();
            const Eigen::Matrix<double,6,6> contourWeight =
                (Eigen::Vector<double,6>() << 1000, 1000, 1000, 40, 40, 40)
                    .finished().asDiagonal();
            const auto objective = [&](const RmpccStateVector &x)
            {
                const auto residual = rmpcc_associated_residuals(
                    x, geometry, 1e-12, arc, transform, derivative);
                return (residual.contour.transpose()
                        * contourWeight * residual.contour)(0)
                    + arc.q_lag() * residual.scalarPoseArcLag
                        * residual.scalarPoseArcLag;
            };
            constexpr double h = 2e-6;
            const double finiteDifference =
                (objective(state + h * direction)
                 - objective(state - h * direction)) / (2.0 * h);
            const double analytic =
                (2.0 * linearized.residual.contour.transpose()
                 * contourWeight * linearized.contourJacobian * direction)(0)
                + 2.0 * arc.q_lag()
                    * linearized.residual.scalarPoseArcLag
                    * (linearized.scalarPoseArcLagJacobian * direction)(0);
            maximumDirectionalDerivativeError = std::max(
                maximumDirectionalDerivativeError,
                std::abs(analytic - finiteDifference)
                / std::max({1.0, std::abs(analytic), std::abs(finiteDifference)}));
        }
        {
            const auto geometry =
                RmpccContourResidualGeometry::ExactAssociatedUnifiedSE3;
            const auto linearized = rmpcc_linearize_associated_residuals(
                state, geometry, 1e-12, 1e-6, arc, transform, derivative);
            Eigen::Vector<double,7> direction;
            direction << -0.15, 0.25, 0.2, -0.35, 0.1, 0.3, -0.2;
            direction.normalize();
            const Eigen::Matrix<double,6,6> contourWeight =
                (Eigen::Vector<double,6>() << 1000, 1000, 1000, 40, 40, 40)
                    .finished().asDiagonal();
            const auto objective = [&](const RmpccStateVector &x)
            {
                const auto residual = rmpcc_associated_residuals(
                    x, geometry, 1e-12, arc, transform, derivative);
                return (residual.contour.transpose()
                        * contourWeight * residual.contour)(0)
                    + (residual.vectorLag.transpose()
                       * lag_weight() * residual.vectorLag)(0);
            };
            constexpr double h = 2e-6;
            const double finiteDifference =
                (objective(state + h * direction)
                 - objective(state - h * direction)) / (2.0 * h);
            const double analytic =
                (2.0 * linearized.residual.contour.transpose()
                 * contourWeight * linearized.contourJacobian * direction)(0)
                + (2.0 * linearized.residual.vectorLag.transpose()
                   * lag_weight() * linearized.vectorLagJacobian * direction)(0);
            maximumDirectionalDerivativeError = std::max(
                maximumDirectionalDerivativeError,
                std::abs(analytic - finiteDifference)
                / std::max({1.0, std::abs(analytic), std::abs(finiteDifference)}));
        }
        require(maximumDirectionalDerivativeError < 1e-6,
                "full objective directional derivative failed");

        std::cout << "arc_samples=" << arc.sample_count()
                  << " arc_interpolation=monotone_linear"
                  << " arc_min_density=" << arc.minimum_density()
                  << " arc_max_density=" << arc.maximum_density()
                  << " arc_total_length=" << arc.total_arc_length()
                  << " arc_max_difference=" << maximumArcDifference
                  << " phase_order=" << phaseOrder
                  << " contour_order=" << contourOrder
                  << " local_exact_order=" << localExactOrder
                  << " arc_order=" << arcOrder
                  << " unified_decoupled_order=" << decoupledOrder
                  << " directional_derivative_error="
                  << maximumDirectionalDerivativeError << '\n';
        return 0;
    }
    catch(const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
