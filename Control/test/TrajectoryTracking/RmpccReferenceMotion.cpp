#include <Control/TrajectoryTracking/SerialLinkRMPCC.h>
#include "RmpccReferenceMotion.h"

#include <Math/MathFunctions.h>

#include <Eigen/Core>

#include <array>
#include <cmath>
#include <iostream>

int main()
{
    using RobotLibrary::Control::RmpccInputVector;
    using RobotLibrary::Control::RmpccParameters;
    using RobotLibrary::Control::RmpccReferenceMotion;
    using RobotLibrary::Control::RmpccStateVector;
    using RobotLibrary::Control::rmpcc_centred_geometric_tangent;
    using RobotLibrary::Control::rmpcc_exact_state_step;
    using RobotLibrary::Control::rmpcc_linearize_stage_consistent_path_velocity_residual;
    using RobotLibrary::Control::rmpcc_stage_consistent_feedforward;
    using RobotLibrary::Control::rmpcc_stage_consistent_path_velocity_residual;
    using RobotLibrary::Math::se3_exponential;
    using RobotLibrary::Math::se3_inverse;
    using RobotLibrary::Math::se3_logarithm;
    using Vector6 = Eigen::Vector<double,6>;

    if(RmpccParameters{}.referenceMotion
       != RmpccReferenceMotion::LegacyTangentProduct)
    {
        std::cerr << "StageConsistent reference motion became the production default.\n";
        return 1;
    }

    Vector6 firstGenerator;
    firstGenerator << 0.12, -0.05, 0.08, 0.35, -0.22, 0.17;
    Vector6 secondGenerator;
    secondGenerator << -0.04, 0.09, 0.03, -0.18, 0.31, 0.11;
    const auto reference = [&](const double progress) -> Eigen::Matrix4d
    {
        const double s = std::clamp(progress, 0.0, 1.0);
        const Eigen::Matrix4d first = se3_exponential(firstGenerator * s);
        const Eigen::Matrix4d second = se3_exponential(secondGenerator * (s * s));
        return first * second;
    };

    const double progress = 0.37;
    const Vector6 tangentA =
        rmpcc_centred_geometric_tangent(progress, 1e-5, reference);
    const Vector6 tangentB =
        rmpcc_centred_geometric_tangent(progress, 1e-5, reference);
    if((tangentA - tangentB).norm() != 0.0)
    {
        std::cerr << "Geometric tangent is not deterministic.\n";
        return 2;
    }

    const std::array<Vector6,3> errors{{
        Vector6::Zero(),
        (Vector6() << 1e-3, -0.7e-3, 0.4e-3, 0.0, 0.0, 0.0).finished(),
        (Vector6() << 0.5e-3, 0.2e-3, -0.8e-3, 5e-3, -3e-3, 4e-3).finished()}};
    double maximumInvariantError = 0.0;
    for(const double dt : {0.001, 0.002, 0.004})
    {
        for(const double rate : {0.02, 0.05, 0.09})
        {
            for(const Vector6 &error : errors)
            {
                RmpccStateVector state = RmpccStateVector::Zero();
                state.head<6>() = error;
                state(6) = progress;
                RmpccInputVector input = RmpccInputVector::Zero();
                input.head<6>() = rmpcc_stage_consistent_feedforward(
                    error, progress, rate, dt, reference);
                input(6) = rate;
                const RmpccStateVector next =
                    rmpcc_exact_state_step(state, input, dt, reference);
                const Vector6 invariant = se3_logarithm(
                    se3_inverse(se3_exponential(error))
                    * se3_exponential(next.head<6>()));
                maximumInvariantError = std::max(maximumInvariantError, invariant.norm());
            }
        }
    }
    if(maximumInvariantError > 1e-11)
    {
        std::cerr << "Stage-consistent feedforward does not preserve E: "
                  << maximumInvariantError << '\n';
        return 3;
    }

    RmpccStateVector state;
    state << 0.012, -0.008, 0.006, 0.031, -0.024, 0.017, progress;
    RmpccInputVector input;
    input << 0.018, -0.011, 0.009, 0.041, -0.027, 0.019, 0.047;
    constexpr double dt = 0.002;
    const auto linearization =
        rmpcc_linearize_stage_consistent_path_velocity_residual(
            state, input, dt, 1e-6, reference);
    if((linearization.inputJacobian.leftCols<6>()
        - Eigen::Matrix<double,6,6>::Identity()).norm() > 1e-15)
    {
        std::cerr << "Stage-consistent residual J_u is not exact identity.\n";
        return 4;
    }

    constexpr double validationStep = 3e-7;
    Eigen::Matrix<double,6,7> expectedState = Eigen::Matrix<double,6,7>::Zero();
    for(int column = 0; column < 7; ++column)
    {
        RmpccStateVector plus = state;
        RmpccStateVector minus = state;
        plus(column) += validationStep;
        minus(column) -= validationStep;
        const double denominator = plus(column) - minus(column);
        expectedState.col(column) =
            (rmpcc_stage_consistent_path_velocity_residual(
                 plus, input, dt, reference)
             - rmpcc_stage_consistent_path_velocity_residual(
                 minus, input, dt, reference)) / denominator;
    }
    RmpccInputVector plusRate = input;
    RmpccInputVector minusRate = input;
    plusRate(6) += validationStep;
    minusRate(6) -= validationStep;
    const Vector6 expectedRate =
        (rmpcc_stage_consistent_path_velocity_residual(
             state, plusRate, dt, reference)
         - rmpcc_stage_consistent_path_velocity_residual(
             state, minusRate, dt, reference)) / (2.0 * validationStep);
    const double stateJacobianError =
        (linearization.stateJacobian - expectedState).norm();
    const double rateJacobianError =
        (linearization.inputJacobian.col(6) - expectedRate).norm();
    if(stateJacobianError > 2e-6 or rateJacobianError > 2e-6)
    {
        std::cerr << "Complete path-velocity residual Jacobian mismatch: state="
                  << stateJacobianError << " rate=" << rateJacobianError << '\n';
        return 5;
    }

    return 0;
}
