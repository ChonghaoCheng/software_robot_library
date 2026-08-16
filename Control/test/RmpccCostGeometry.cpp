#include <Control/RmpccCostGeometry.h>
#include <Control/RmpccResidualLinearization.h>

#include <Eigen/Core>

#include <iostream>

int main()
{
    using RobotLibrary::Control::RmpccLagGeometry;
    using RobotLibrary::Control::rmpcc_component_scaled_weight;
    using RobotLibrary::Control::rmpcc_error_projection;
    using RobotLibrary::Control::rmpcc_decoupled_cost_weight;
    using RobotLibrary::Control::rmpcc_decoupled_error;
    using RobotLibrary::Control::rmpcc_decoupled_error_projection;
    using RobotLibrary::Control::rmpcc_decoupled_residuals;

    Eigen::Matrix<double,6,1> tangent;
    tangent << 0.20, -0.08, 0.04, 0.35, 0.12, -0.18;
    Eigen::Matrix<double,6,6> metric = Eigen::Matrix<double,6,6>::Identity();
    metric.diagonal().tail<3>().setConstant(0.04);
    constexpr double epsilon = 1e-12;

    const auto full = rmpcc_error_projection(
        tangent, metric, RmpccLagGeometry::FullScrew, epsilon);
    if(full.lag.block<3,3>(0,3).norm() < 1e-6
       or full.lag.block<3,3>(3,0).norm() < 1e-6)
    {
        std::cerr << "Full screw lag unexpectedly lacks translation/rotation coupling.\n";
        return 1;
    }

    const auto split = rmpcc_error_projection(
        tangent, metric, RmpccLagGeometry::SplitTranslationRotation, epsilon);
    if(split.lag.block<3,3>(0,3).norm() > 1e-14
       or split.lag.block<3,3>(3,0).norm() > 1e-14)
    {
        std::cerr << "Split lag retains cross-subspace coupling.\n";
        return 2;
    }

    const auto translation = rmpcc_error_projection(
        tangent, metric, RmpccLagGeometry::TranslationOnly, epsilon);
    if(translation.lag.bottomRows<3>().norm() > 1e-14
       or (translation.contour.block<3,3>(3,3) - Eigen::Matrix3d::Identity()).norm() > 1e-14)
    {
        std::cerr << "Translation-only lag does not leave rotation as contour error.\n";
        return 3;
    }

    const auto rotation = rmpcc_error_projection(
        tangent, metric, RmpccLagGeometry::RotationOnly, epsilon);
    if(rotation.lag.topRows<3>().norm() > 1e-14
       or (rotation.contour.block<3,3>(0,0) - Eigen::Matrix3d::Identity()).norm() > 1e-14)
    {
        std::cerr << "Rotation-only lag does not leave translation as contour error.\n";
        return 4;
    }

    if((split.lag * split.lag - split.lag).norm() > 1e-12
       or (split.lag + split.contour
           - Eigen::Matrix<double,6,6>::Identity()).norm() > 1e-14)
    {
        std::cerr << "Split contour/lag projectors are inconsistent.\n";
        return 5;
    }

    Eigen::Matrix<double,6,1> translationOnlyTangent = tangent;
    translationOnlyTangent.tail<3>().setZero();
    const auto splitWithoutRotation = rmpcc_error_projection(
        translationOnlyTangent, metric,
        RmpccLagGeometry::SplitTranslationRotation, epsilon);
    const auto translationWithoutRotation = rmpcc_error_projection(
        translationOnlyTangent, metric,
        RmpccLagGeometry::TranslationOnly, epsilon);
    if((splitWithoutRotation.lag - translationWithoutRotation.lag).norm() > 1e-14)
    {
        std::cerr << "A zero angular tangent should make split and translation-only lag equal.\n";
        return 6;
    }

    Eigen::Matrix<double,6,6> weight = Eigen::Matrix<double,6,6>::Identity();
    weight.diagonal() << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
    const auto translationCostOnly =
        rmpcc_component_scaled_weight(weight, 1.0, 0.0);
    if((translationCostOnly.diagonal().head<3>() - weight.diagonal().head<3>()).norm() > 1e-14
       or translationCostOnly.bottomRows<3>().norm() > 1e-14
       or translationCostOnly.rightCols<3>().norm() > 1e-14)
    {
        std::cerr << "Component scaling failed to remove only the rotational cost.\n";
        return 7;
    }

    const auto rotationCostOnly =
        rmpcc_component_scaled_weight(weight, 0.0, 1.0);
    if(rotationCostOnly.topRows<3>().norm() > 1e-14
       or rotationCostOnly.leftCols<3>().norm() > 1e-14
       or (rotationCostOnly.diagonal().tail<3>() - weight.diagonal().tail<3>()).norm() > 1e-14)
    {
        std::cerr << "Component scaling failed to remove only the translational cost.\n";
        return 8;
    }

    Eigen::Matrix<double,6,1> se3Error;
    se3Error << 0.03, -0.02, 0.01, 0.4, -0.2, 0.15;
    const Eigen::Matrix<double,6,1> decoupled = rmpcc_decoupled_error(se3Error);
    const Eigen::Matrix4d relative =
        RobotLibrary::Math::se3_exponential(se3Error);
    if((decoupled.head<3>() - relative.block<3,1>(0,3)).norm() > 1e-14
       or (decoupled.tail<3>() - se3Error.tail<3>()).norm() > 1e-12)
    {
        std::cerr << "Decoupled residual is not Cartesian position plus SO(3) error.\n";
        return 9;
    }

    const auto decoupledWeight = rmpcc_decoupled_cost_weight(
        tangent, weight, 2.0 * weight, epsilon);
    Eigen::Matrix<double,6,1> changedAngularTangent = tangent;
    changedAngularTangent.tail<3>() *= -7.0;
    const auto changedAngularWeight = rmpcc_decoupled_cost_weight(
        changedAngularTangent, weight, 2.0 * weight, epsilon);
    if((decoupledWeight - changedAngularWeight).norm() > 1e-14)
    {
        std::cerr << "Independent SO(3) cost contains a hidden angular phase projector.\n";
        return 10;
    }

    Eigen::Matrix<double,6,1> pureRotationTangent = Eigen::Matrix<double,6,1>::Zero();
    pureRotationTangent.tail<3>() << 0.7, -0.4, 0.3;
    const auto pureRotationProjection = rmpcc_decoupled_error_projection(
        pureRotationTangent, epsilon);
    if(pureRotationProjection.lag.norm() > 1e-14
       or (pureRotationProjection.contour
           - Eigen::Matrix<double,6,6>::Identity()).norm() > 1e-14)
    {
        std::cerr << "Pure rotation must have zero translational lag projector.\n";
        return 11;
    }
    Eigen::Matrix<double,6,1> arbitraryError;
    arbitraryError << 0.13, -0.09, 0.07, 0.5, -0.3, 0.2;
    if((pureRotationProjection.lag * arbitraryError).norm() > 1e-14
       or ((pureRotationProjection.contour * arbitraryError).head<3>()
           - arbitraryError.head<3>()).norm() > 1e-14)
    {
        std::cerr << "Pure-rotation geometry classified Cartesian error as lag.\n";
        return 12;
    }

    Eigen::Matrix<double,6,1> smallTangent = pureRotationTangent;
    smallTangent.head<3>() << 0.0, 1e-8, 0.0;
    const double smallLagNorm =
        rmpcc_decoupled_error_projection(smallTangent, 1e-9).lag.norm();
    smallTangent.head<3>() *= 0.5;
    const double smallerLagNorm =
        rmpcc_decoupled_error_projection(smallTangent, 1e-9).lag.norm();
    if(not (smallerLagNorm < smallLagNorm and smallLagNorm < 1e-6)
       or rmpcc_decoupled_error_projection(smallTangent, 1e-9).lag(0,0) != 0.0)
    {
        std::cerr << "Near-zero tangent does not approach zero lag continuously.\n";
        return 13;
    }

    RobotLibrary::Control::RmpccStateVector pureRotationState;
    pureRotationState << arbitraryError, 0.37;
    const auto pureRotationResidual = rmpcc_decoupled_residuals(
        pureRotationState, epsilon,
        [&](const double) { return pureRotationTangent; });
    const Eigen::Matrix<double,6,1> mappedError =
        rmpcc_decoupled_error(pureRotationState.head<6>());
    if((pureRotationResidual.contour
        - pureRotationProjection.contour * mappedError).norm() > 1e-14
       or (pureRotationResidual.lag
           - pureRotationProjection.lag * mappedError).norm() > 1e-14)
    {
        std::cerr << "Optimizer residual and diagnostic projection disagree.\n";
        return 14;
    }
    const auto pureRotationWeight = rmpcc_decoupled_cost_weight(
        pureRotationTangent, weight, 2.0 * weight, epsilon);
    const double weightedError =
        (mappedError.transpose() * pureRotationWeight * mappedError)(0);
    const double residualCost =
        (pureRotationResidual.contour.transpose()
         * weight * pureRotationResidual.contour)(0)
        + (pureRotationResidual.lag.transpose()
           * (2.0 * weight) * pureRotationResidual.lag)(0);
    if(std::abs(weightedError - residualCost) > 1e-12)
    {
        std::cerr << "Frozen optimizer weight and residual diagnostics disagree.\n";
        return 15;
    }

    return 0;
}
