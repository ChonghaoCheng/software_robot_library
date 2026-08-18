/**
 * @file TrajectoryTracking/test/RmpccResidualNonExpansion.cpp
 * @brief F2-01/F2-02 phase-association conditioning robustness gate.
 */

#include "RmpccPhaseResidual.h"

#include <Math/MathFunctions.h>

#include <Eigen/Core>
#include <Eigen/SVD>

#include <array>
#include <cmath>
#include <fstream>
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
constexpr double kLinearTangent = 0.176;
constexpr double kFiniteDifferenceStep = 1e-7;
constexpr double kPureProgressEpsilon = 1e-3;
constexpr double kRotationOnlyError = 1e-3;

double gAnisotropy = 1.0;
double gCharacteristicLength = 0.20;
Eigen::Vector3d gRotationAxis = Eigen::Vector3d::UnitZ();

struct Axis
{
    const char *label;
    Eigen::Vector3d direction;
};

std::vector<Axis> axes()
{
    // The R0--R7 directions are fixed pseudo-random draws stored explicitly,
    // so the robustness gate is deterministic across standard libraries.
    std::vector<Axis> result{
        {"X", Eigen::Vector3d::UnitX()},
        {"Y", Eigen::Vector3d::UnitY()},
        {"Z", Eigen::Vector3d::UnitZ()},
        {"R0", Eigen::Vector3d(1.0, 2.0, 3.0)},
        {"R1", Eigen::Vector3d(-2.0, 1.0, 4.0)},
        {"R2", Eigen::Vector3d(3.0, -4.0, 1.0)},
        {"R3", Eigen::Vector3d(-1.0, -3.0, 2.0)},
        {"R4", Eigen::Vector3d(4.0, 2.0, -3.0)},
        {"R5", Eigen::Vector3d(2.0, -5.0, -1.0)},
        {"R6", Eigen::Vector3d(-3.0, 4.0, 5.0)},
        {"R7", Eigen::Vector3d(5.0, -1.0, -2.0)}};
    for(Axis &axis : result) axis.direction.normalize();
    return result;
}

Eigen::Vector<double,6> screw()
{
    Eigen::Vector<double,6> value;
    value.head<3>() = kLinearTangent * Eigen::Vector3d::UnitX();
    value.tail<3>() =
        gAnisotropy * kLinearTangent * gRotationAxis;
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
        gCharacteristicLength * gCharacteristicLength);
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
        gCharacteristicLength, reference, tangent);
}

RobotLibrary::Control::RmpccPhaseResidualLinearization linearization(
    const RmpccStateVector &state,
    const RmpccPhaseAssociation association)
{
    return RobotLibrary::Control::rmpcc_linearize_phase_residuals(
        state, metric(), association, kRegularization, kTolerance,
        gCharacteristicLength, kFiniteDifferenceStep, reference, tangent);
}

double spectral_norm(const Eigen::MatrixXd &matrix)
{
    return Eigen::JacobiSVD<Eigen::MatrixXd>(matrix).singularValues()(0);
}

struct Row
{
    double ratio = 0.0;
    std::string axis;
    Eigen::Vector3d axisDirection = Eigen::Vector3d::Zero();
    double rotationLength = 0.0;
    RmpccPhaseAssociation association = RmpccPhaseAssociation::MetricScrew;
    double contourProjectionNorm = 0.0;
    double rotationFromPositionNorm = 0.0;
    double pureProgressError = 0.0;
    double pureProgressContour = 0.0;
    double pureRotationPhaseShift = 0.0;
    double phaseDenominator = 0.0;
};

Row measure(
    const double ratio,
    const Axis &axis,
    const double rotationLength,
    const RmpccPhaseAssociation association)
{
    gAnisotropy = ratio;
    gRotationAxis = axis.direction;
    gCharacteristicLength = rotationLength;

    RmpccStateVector onPath = RmpccStateVector::Zero();
    onPath(6) = 0.5;
    const auto local = linearization(onPath, association);

    RmpccStateVector progressed = RmpccStateVector::Zero();
    progressed.head<6>() = screw() * kPureProgressEpsilon;
    progressed(6) = 0.5;
    const auto progressResidual = residual(progressed, association);

    RmpccStateVector rotated = RmpccStateVector::Zero();
    rotated.segment<3>(3) = kRotationOnlyError * axis.direction;
    rotated(6) = 0.5;
    const auto rotationResidual = residual(rotated, association);

    Row row;
    row.ratio = ratio;
    row.axis = axis.label;
    row.axisDirection = axis.direction;
    row.rotationLength = rotationLength;
    row.association = association;
    row.contourProjectionNorm = spectral_norm(
        local.contourJacobian.leftCols<6>());
    row.rotationFromPositionNorm = spectral_norm(
        local.contourJacobian.block<3,3>(3,0));
    row.pureProgressError = std::abs(
        progressResidual.phaseCorrection - kPureProgressEpsilon);
    row.pureProgressContour = progressResidual.contour.norm();
    row.pureRotationPhaseShift = std::abs(rotationResidual.phaseCorrection);
    row.phaseDenominator = local.residual.phaseDenominator;
    return row;
}

bool pure_progress_gate(
    const double ratio,
    const Axis &axis,
    const double rotationLength)
{
    gAnisotropy = ratio;
    gRotationAxis = axis.direction;
    gCharacteristicLength = rotationLength;
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
                      << ratio << ", axis=" << axis.label
                      << ", l_R=" << rotationLength
                      << ", epsilon=" << epsilon
                      << ": phase_error=" << phaseError
                      << ", contour=" << value.contour.norm()
                      << ", O(epsilon^2) bound=" << bound << "\n";
            ok = false;
        }
    }
    return ok;
}

struct LengthSummary
{
    double rotationLength = 0.0;
    double maximumProjectionNorm = 0.0;
    double maximumCoupling = 0.0;
    double maximumProgressError = 0.0;
    double maximumProgressContour = 0.0;
    double maximumRotationPhaseShift = 0.0;
    double minimumDenominator = std::numeric_limits<double>::infinity();
};

void write_csv(std::ostream &output, const std::vector<Row> &rows)
{
    output << std::setprecision(17);
    output << "ratio,axis,axis_x,axis_y,axis_z,l_R,association,Pc_norm,"
              "d_contour_R_d_e_p_norm,pure_progress_error,"
              "pure_progress_contour,pure_rotation_abs_delta_s,"
              "phase_denominator\n";
    for(const Row &row : rows)
    {
        output << row.ratio << ',' << row.axis << ','
               << row.axisDirection.x() << ',' << row.axisDirection.y() << ','
               << row.axisDirection.z() << ',' << row.rotationLength << ','
               << name(row.association) << ','
               << row.contourProjectionNorm << ','
               << row.rotationFromPositionNorm << ','
               << row.pureProgressError << ','
               << row.pureProgressContour << ','
               << row.pureRotationPhaseShift << ','
               << row.phaseDenominator << '\n';
    }
}

} // namespace

int main(int argc, char **argv)
{
    if(argc > 2)
    {
        std::cerr << "Usage: " << argv[0] << " [output.csv]\n";
        return 2;
    }
    constexpr std::array<double,5> ratios{{1.0, 5.0, 10.0, 20.0, 50.0}};
    constexpr std::array<double,4> rotationLengths{{0.05, 0.10, 0.20, 0.40}};
    constexpr std::array<RmpccPhaseAssociation,3> associations{{
        RmpccPhaseAssociation::MetricScrew,
        RmpccPhaseAssociation::TaskPointXYZ,
        RmpccPhaseAssociation::TaskPoseFeature}};
    const std::vector<Axis> rotationAxes = axes();

    std::vector<Row> rows;
    std::vector<LengthSummary> summaries;
    int failures = 0;
    for(const double rotationLength : rotationLengths)
    {
        LengthSummary summary;
        summary.rotationLength = rotationLength;
        for(const Axis &axis : rotationAxes)
        {
            const double axisFeatureFactor =
                1.0 + axis.direction.z() * axis.direction.z();
            const double couplingBound =
                1.0 / (2.0 * rotationLength * std::sqrt(axisFeatureFactor));
            const double featureJacobianNorm =
                std::max(1.0, rotationLength * std::sqrt(2.0));
            const double projectionBound = featureJacobianNorm
                * std::max(1.0,
                           1.0 / (rotationLength * std::sqrt(axisFeatureFactor)));
            for(const double ratio : ratios)
            {
                failures += pure_progress_gate(
                    ratio, axis, rotationLength) ? 0 : 1;
                for(const auto association : associations)
                {
                    const Row row = measure(
                        ratio, axis, rotationLength, association);
                    rows.push_back(row);
                    if(association == RmpccPhaseAssociation::TaskPointXYZ)
                    {
                        const double relativeError = std::abs(
                            row.rotationFromPositionNorm - ratio) / ratio;
                        if(relativeError > 1e-6)
                        {
                            std::cerr << "TaskPointXYZ known gain mismatch at ratio="
                                      << ratio << ", axis=" << axis.label
                                      << ": measured="
                                      << row.rotationFromPositionNorm
                                      << ", expected=" << ratio << "\n";
                            ++failures;
                        }
                    }
                    if(association == RmpccPhaseAssociation::TaskPoseFeature)
                    {
                        summary.maximumProjectionNorm = std::max(
                            summary.maximumProjectionNorm,
                            row.contourProjectionNorm);
                        summary.maximumCoupling = std::max(
                            summary.maximumCoupling,
                            row.rotationFromPositionNorm);
                        summary.maximumProgressError = std::max(
                            summary.maximumProgressError,
                            row.pureProgressError);
                        summary.maximumProgressContour = std::max(
                            summary.maximumProgressContour,
                            row.pureProgressContour);
                        summary.maximumRotationPhaseShift = std::max(
                            summary.maximumRotationPhaseShift,
                            row.pureRotationPhaseShift);
                        summary.minimumDenominator = std::min(
                            summary.minimumDenominator,
                            row.phaseDenominator);
                        if(row.rotationFromPositionNorm > couplingBound + 1e-5)
                        {
                            std::cerr << "TaskPoseFeature coupling bound failed: ratio="
                                      << ratio << ", axis=" << axis.label
                                      << ", l_R=" << rotationLength
                                      << ", measured="
                                      << row.rotationFromPositionNorm
                                      << ", bound=" << couplingBound << "\n";
                            ++failures;
                        }
                        if(row.contourProjectionNorm > projectionBound + 1e-5)
                        {
                            std::cerr << "TaskPoseFeature projection bound failed: ratio="
                                      << ratio << ", axis=" << axis.label
                                      << ", l_R=" << rotationLength
                                      << ", measured=" << row.contourProjectionNorm
                                      << ", bound=" << projectionBound << "\n";
                            ++failures;
                        }
                        const double rotationShiftBound =
                            2.0 * kRotationOnlyError / kLinearTangent;
                        if(not std::isfinite(row.pureRotationPhaseShift)
                           or row.pureRotationPhaseShift > rotationShiftBound)
                        {
                            std::cerr << "TaskPoseFeature rotation-only phase shift "
                                         "is unbounded: ratio="
                                      << ratio << ", axis=" << axis.label
                                      << ", l_R=" << rotationLength
                                      << ", measured="
                                      << row.pureRotationPhaseShift
                                      << ", bound=" << rotationShiftBound << "\n";
                            ++failures;
                        }
                    }
                }
            }
        }
        summaries.push_back(summary);
    }

    if(argc == 2)
    {
        std::ofstream output(argv[1]);
        if(not output)
        {
            std::cerr << "Cannot open CSV output: " << argv[1] << '\n';
            return 3;
        }
        write_csv(output, rows);
    }

    std::cout << std::setprecision(17);
    std::cout << "l_R,max_TaskPoseFeature_Pc,max_TaskPoseFeature_coupling,"
                 "max_pure_progress_error,max_pure_progress_contour,"
                 "max_pure_rotation_abs_delta_s,min_phase_denominator\n";
    for(const LengthSummary &summary : summaries)
    {
        std::cout << summary.rotationLength << ','
                  << summary.maximumProjectionNorm << ','
                  << summary.maximumCoupling << ','
                  << summary.maximumProgressError << ','
                  << summary.maximumProgressContour << ','
                  << summary.maximumRotationPhaseShift << ','
                  << summary.minimumDenominator << '\n';
    }
    std::cout << "rows," << rows.size() << ",axes," << rotationAxes.size()
              << ",ratios," << ratios.size() << ",rotation_lengths,"
              << rotationLengths.size() << '\n';
    return failures == 0 ? 0 : 1;
}
