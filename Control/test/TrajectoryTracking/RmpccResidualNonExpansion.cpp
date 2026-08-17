/**
 * @file TrajectoryTracking/test/RmpccResidualNonExpansion.cpp
 * @brief F2-01 phase-association conditioning and consistency gate.
 */

#include "RmpccPhaseResidual.h"

#include <Math/MathFunctions.h>

#include <Eigen/Core>
#include <Eigen/SVD>

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using RobotLibrary::Control::RmpccPhaseAssociation;
using RobotLibrary::Control::RmpccStateVector;

constexpr double kRegularization = 1e-10;
constexpr double kTolerance = 1e-12;
constexpr double kCharacteristicLength = 0.20;
constexpr double kLinearTangent = 0.176;
constexpr double kFiniteDifferenceStep = 1e-7;
constexpr double kPureProgressEpsilon = 1e-3;
constexpr double kRotationOnlyError = 1e-3;

double gAnisotropy = 1.0;

Eigen::Vector<double,6> screw()
{
    Eigen::Vector<double,6> value;
    value << kLinearTangent, 0.0, 0.0,
             0.0, 0.0, gAnisotropy * kLinearTangent;
    return value;
}

Eigen::Matrix4d reference(const double progress)
{
    return RobotLibrary::Math::se3_exponential(screw() * progress);
}

Eigen::Vector<double,6> tangent(const double) { return screw(); }

Eigen::Matrix<double,6,6> metric()
{
    Eigen::Matrix<double,6,6> value = Eigen::Matrix<double,6,6>::Identity();
    value.diagonal().tail<3>().setConstant(
        kCharacteristicLength * kCharacteristicLength);
    return value;
}

const char *name(const RmpccPhaseAssociation association)
{
    switch(association)
    {
        case RmpccPhaseAssociation::MetricScrew: return "MetricScrew";
        case RmpccPhaseAssociation::TaskPointXYZ: return "TaskPointXYZ";
        case RmpccPhaseAssociation::TaskPoseFeature: return "TaskPoseFeature";
    }
    return "Unknown";
}

RobotLibrary::Control::RmpccPhaseResiduals residual(
    const RmpccStateVector &state,
    const RmpccPhaseAssociation association)
{
    return RobotLibrary::Control::rmpcc_phase_residuals(
        state, metric(), association, kRegularization, kTolerance,
        kCharacteristicLength, reference, tangent);
}

RobotLibrary::Control::RmpccPhaseResidualLinearization linearization(
    const RmpccStateVector &state,
    const RmpccPhaseAssociation association)
{
    return RobotLibrary::Control::rmpcc_linearize_phase_residuals(
        state, metric(), association, kRegularization, kTolerance,
        kCharacteristicLength, kFiniteDifferenceStep, reference, tangent);
}

double spectral_norm(const Eigen::MatrixXd &matrix)
{
    return Eigen::JacobiSVD<Eigen::MatrixXd>(matrix).singularValues()(0);
}

struct Row
{
    double ratio = 0.0;
    RmpccPhaseAssociation association = RmpccPhaseAssociation::MetricScrew;
    double contourProjectionNorm = 0.0;
    double rotationFromPositionNorm = 0.0;
    double pureProgressError = 0.0;
    double pureProgressContour = 0.0;
    double phaseDenominator = 0.0;
};

Row measure(const double ratio, const RmpccPhaseAssociation association)
{
    gAnisotropy = ratio;
    RmpccStateVector onPath = RmpccStateVector::Zero();
    onPath(6) = 0.5;
    const auto local = linearization(onPath, association);

    RmpccStateVector progressed = RmpccStateVector::Zero();
    progressed.head<6>() = screw() * kPureProgressEpsilon;
    progressed(6) = 0.5;
    const auto progressResidual = residual(progressed, association);

    Row row;
    row.ratio = ratio;
    row.association = association;
    row.contourProjectionNorm = spectral_norm(
        local.contourJacobian.leftCols<6>());
    row.rotationFromPositionNorm = spectral_norm(
        local.contourJacobian.block<3,3>(3,0));
    row.pureProgressError = std::abs(
        progressResidual.phaseCorrection - kPureProgressEpsilon);
    row.pureProgressContour = progressResidual.contour.norm();
    row.phaseDenominator = local.residual.phaseDenominator;
    return row;
}

bool pure_progress_gate(const double ratio)
{
    gAnisotropy = ratio;
    bool ok = true;
    for(const double epsilon : {2e-3, 1e-3, 5e-4})
    {
        RmpccStateVector progressed = RmpccStateVector::Zero();
        progressed.head<6>() = screw() * epsilon;
        progressed(6) = 0.5;
        const auto value = residual(
            progressed, RmpccPhaseAssociation::TaskPoseFeature);
        const double phaseError = std::abs(value.phaseCorrection - epsilon);
        const double bound = 10.0 * epsilon * epsilon;
        if(phaseError > bound or value.contour.norm() > bound)
        {
            std::cerr << "TaskPoseFeature pure-progress consistency failed at ratio="
                      << ratio << ", epsilon=" << epsilon
                      << ": phase_error=" << phaseError
                      << ", contour=" << value.contour.norm()
                      << ", O(epsilon^2) bound=" << bound << "\n";
            ok = false;
        }
    }
    return ok;
}

bool rotation_only_gate(const double ratio)
{
    gAnisotropy = ratio;
    RmpccStateVector state = RmpccStateVector::Zero();
    state(5) = kRotationOnlyError;
    state(6) = 0.5;
    const auto value = residual(state, RmpccPhaseAssociation::TaskPoseFeature);
    const double conservativeBound =
        2.0 * kRotationOnlyError / kLinearTangent;
    if(not std::isfinite(value.phaseCorrection)
       or std::abs(value.phaseCorrection) > conservativeBound)
    {
        std::cerr << "TaskPoseFeature rotation-only phase shift is unbounded at ratio="
                  << ratio << ": delta_s=" << value.phaseCorrection
                  << ", bound=" << conservativeBound << "\n";
        return false;
    }
    return true;
}

} // namespace

int main()
{
    constexpr std::array<double,6> ratios{{1.0, 2.0, 5.0, 10.0, 20.0, 50.0}};
    constexpr std::array<RmpccPhaseAssociation,3> associations{{
        RmpccPhaseAssociation::MetricScrew,
        RmpccPhaseAssociation::TaskPointXYZ,
        RmpccPhaseAssociation::TaskPoseFeature}};

    std::vector<Row> rows;
    int failures = 0;
    double maximumFeatureCoupling = 0.0;
    for(const double ratio : ratios)
    {
        failures += pure_progress_gate(ratio) ? 0 : 1;
        failures += rotation_only_gate(ratio) ? 0 : 1;
        for(const auto association : associations)
        {
            const Row row = measure(ratio, association);
            rows.push_back(row);
            if(association == RmpccPhaseAssociation::TaskPointXYZ)
            {
                const double relativeError =
                    std::abs(row.rotationFromPositionNorm - ratio) / ratio;
                if(relativeError > 1e-6)
                {
                    std::cerr << "TaskPointXYZ known gain mismatch at ratio="
                              << ratio << ": measured="
                              << row.rotationFromPositionNorm
                              << ", expected=" << ratio << "\n";
                    ++failures;
                }
            }
            if(association == RmpccPhaseAssociation::TaskPoseFeature)
            {
                maximumFeatureCoupling = std::max(
                    maximumFeatureCoupling, row.rotationFromPositionNorm);
            }
        }
    }

    const double featureCouplingBound =
        1.0 / (2.0 * std::sqrt(2.0) * kCharacteristicLength);
    if(maximumFeatureCoupling > featureCouplingBound + 1e-5)
    {
        std::cerr << "TaskPoseFeature position-to-rotation coupling exceeds "
                  << "the analytic bound: measured=" << maximumFeatureCoupling
                  << ", bound=" << featureCouplingBound << "\n";
        ++failures;
    }

    std::cout << std::setprecision(17);
    std::cout << "ratio,association,Pc_norm,d_contour_R_d_e_p_norm,"
                 "pure_progress_error,pure_progress_contour,phase_denominator\n";
    for(const Row &row : rows)
    {
        std::cout << row.ratio << ',' << name(row.association) << ','
                  << row.contourProjectionNorm << ','
                  << row.rotationFromPositionNorm << ','
                  << row.pureProgressError << ','
                  << row.pureProgressContour << ','
                  << row.phaseDenominator << '\n';
    }
    std::cout << "TaskPoseFeature coupling bound," << featureCouplingBound
              << ", measured maximum," << maximumFeatureCoupling << '\n';
    return failures == 0 ? 0 : 1;
}
