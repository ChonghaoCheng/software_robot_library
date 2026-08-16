#include "RmpccPrediction.h"
#include <Math/MathFunctions.h>

#include <Eigen/Core>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

namespace {

template<typename StepFunction, typename LinearizeFunction>
bool progress_rate_boundary_jacobian_matches(
    const char *name,
    const RobotLibrary::Control::RmpccStateVector &baseState,
    const RobotLibrary::Control::RmpccInputVector &baseInput,
    const double dt,
    StepFunction step,
    LinearizeFunction linearize)
{
    using RobotLibrary::Control::RmpccInputVector;
    using RobotLibrary::Control::RmpccStateVector;
    constexpr double h = 1e-6;
    for(const double progress :
        std::array<double,5>{0.0, 0.5 * dt * h, 0.4, 1.0 - 0.5 * dt * h, 1.0})
    {
        RmpccStateVector state = baseState;
        state(6) = progress;
        RmpccInputVector input = baseInput;
        input(6) = 0.0;

        const double minimumRate = -progress / dt;
        const double maximumRate = (1.0 - progress) / dt;
        RmpccInputVector minus = input;
        RmpccInputVector plus = input;
        minus(6) = std::clamp(input(6) - h, minimumRate, maximumRate);
        plus(6) = std::clamp(input(6) + h, minimumRate, maximumRate);
        const double denominator = plus(6) - minus(6);
        const RmpccStateVector expected =
            (step(state, plus) - step(state, minus)) / denominator;
        const RmpccStateVector actual = linearize(state, input).inputJacobian.col(6);
        const double error = (actual - expected).norm();
        if(error > 1e-7)
        {
            std::cerr << name << " progress-rate Jacobian mismatch at s="
                      << progress << ": " << error << '\n';
            return false;
        }
    }
    return true;
}

} // namespace

int main()
{
    using RobotLibrary::Control::RmpccInputVector;
    using RobotLibrary::Control::RmpccStateVector;
    using RobotLibrary::Control::rmpcc_exact_state_step;
    using RobotLibrary::Control::rmpcc_additive_state_step;
    using RobotLibrary::Control::rmpcc_error_coordinate_path_tangent;
    using RobotLibrary::Control::rmpcc_linearize_additive_state_step;
    using RobotLibrary::Control::rmpcc_linearize_exact_state_step;
    using RobotLibrary::Control::rmpcc_transport_reference_tangent;
    using RobotLibrary::Math::se3_exponential;
    using RobotLibrary::Math::se3_inverse;
    using RobotLibrary::Math::se3_logarithm;

    Eigen::Matrix<double,6,1> referenceGenerator;
    referenceGenerator << 0.16, -0.04, 0.08, 0.22, -0.18, 0.11;
    const auto reference = [&referenceGenerator](const double progress)
    {
        return se3_exponential(progress * referenceGenerator);
    };

    RmpccStateVector state;
    state << 0.05, -0.03, 0.02, 0.18, -0.12, 0.09, 0.31;
    RmpccInputVector input;
    input << 0.08, -0.02, 0.04, -0.03, 0.06, 0.02, 0.27;
    const double dt = 0.012;

    const RmpccStateVector next = rmpcc_exact_state_step(state, input, dt, reference);
    const double nextProgress = state(6) + dt * input(6);
    const Eigen::Matrix4d expectedRelative =
        se3_inverse(reference(nextProgress))
        * reference(state(6))
        * se3_exponential(state.head<6>())
        * se3_exponential(dt * input.head<6>());
    if((next.head<6>() - se3_logarithm(expectedRelative)).norm() > 1e-12
       or std::abs(next(6) - nextProgress) > 1e-14)
    {
        std::cerr << "Exact RMPCC group rollout does not match the discrete definition.\n";
        return 1;
    }

    const auto linearization =
        rmpcc_linearize_exact_state_step(state, input, dt, 1e-6, reference);
    RmpccStateVector statePerturbation;
    statePerturbation << 2e-6, -1e-6, 1.5e-6, -1.2e-6, 0.8e-6, 1.1e-6, -1.0e-6;
    RmpccInputVector inputPerturbation;
    inputPerturbation << -1.1e-6, 0.7e-6, -0.9e-6, 1.3e-6, -0.8e-6, 0.6e-6, 0.9e-6;
    const RmpccStateVector nonlinear =
        rmpcc_exact_state_step(state + statePerturbation,
                               input + inputPerturbation,
                               dt,
                               reference);
    const RmpccStateVector predicted =
        linearization.nominalNext
        + linearization.stateJacobian * statePerturbation
        + linearization.inputJacobian * inputPerturbation;
    const double residual = (nonlinear - predicted).norm();
    if(residual > 2e-10)
    {
        std::cerr << "RMPCC RTI linearization residual: " << residual << '\n';
        return 2;
    }

    // The contour/lag direction is the derivative of logarithmic error
    // coordinates, not the untransported reference body tangent.
    const double tangentStep = 1e-7;
    const double tangentProgress = state(6);
    const Eigen::Matrix4d fixedEndpoint =
        reference(tangentProgress) * se3_exponential(state.head<6>());
    const Eigen::Matrix<double,6,1> errorPlus = se3_logarithm(
        se3_inverse(reference(tangentProgress + tangentStep)) * fixedEndpoint);
    const Eigen::Matrix<double,6,1> errorMinus = se3_logarithm(
        se3_inverse(reference(tangentProgress - tangentStep)) * fixedEndpoint);
    const Eigen::Matrix<double,6,1> finiteDifferencePathDirection =
        -(errorPlus - errorMinus) / (2.0 * tangentStep);
    const Eigen::Matrix<double,6,1> errorCoordinateTangent =
        rmpcc_error_coordinate_path_tangent(state.head<6>(), referenceGenerator);
    if((finiteDifferencePathDirection - errorCoordinateTangent).norm() > 2e-7)
    {
        std::cerr << "RMPCC error-coordinate tangent mismatch: "
                  << (finiteDifferencePathDirection - errorCoordinateTangent).norm()
                  << '\n';
        return 3;
    }

    // The transported path velocity must preserve a non-zero relative pose
    // under the exact group rollout.
    RmpccInputVector matchedInput = RmpccInputVector::Zero();
    matchedInput(6) = input(6);
    matchedInput.head<6>() = rmpcc_transport_reference_tangent(
        state.head<6>(), referenceGenerator) * matchedInput(6);
    const RmpccStateVector matchedNext =
        rmpcc_exact_state_step(state, matchedInput, dt, reference);
    if((matchedNext.head<6>() - state.head<6>()).norm() > 2e-12)
    {
        std::cerr << "Transported reference tangent does not preserve relative pose: "
                  << (matchedNext.head<6>() - state.head<6>()).norm() << '\n';
        return 4;
    }

    const auto referenceTangent = [&referenceGenerator](double)
    {
        return referenceGenerator;
    };
    const RmpccStateVector additiveMatchedNext =
        rmpcc_additive_state_step(state, matchedInput, dt, referenceTangent);
    if((additiveMatchedNext.head<6>() - state.head<6>()).norm() > 2e-12)
    {
        std::cerr << "Matched path velocity does not preserve additive error: "
                  << (additiveMatchedNext.head<6>() - state.head<6>()).norm() << '\n';
        return 6;
    }

    const auto exactStep = [&](const RmpccStateVector &x, const RmpccInputVector &u)
    {
        return rmpcc_exact_state_step(x, u, dt, reference);
    };
    const auto exactLinearize = [&](const RmpccStateVector &x, const RmpccInputVector &u)
    {
        return rmpcc_linearize_exact_state_step(x, u, dt, 1e-6, reference);
    };
    const auto additiveStep = [&](const RmpccStateVector &x, const RmpccInputVector &u)
    {
        return rmpcc_additive_state_step(x, u, dt, referenceTangent);
    };
    const auto additiveLinearize = [&](const RmpccStateVector &x, const RmpccInputVector &u)
    {
        return rmpcc_linearize_additive_state_step(
            x, u, dt, 1e-6, referenceTangent);
    };
    if(not progress_rate_boundary_jacobian_matches(
           "ExactSE3", state, input, dt, exactStep, exactLinearize)
       or not progress_rate_boundary_jacobian_matches(
           "Additive", state, input, dt, additiveStep, additiveLinearize))
    {
        return 8;
    }

    // Matrix lag weighting is backward compatible with the former scalar
    // expression when the metric and lag weights are isotropic.
    const Eigen::Matrix<double,6,1> lagTangent = errorCoordinateTangent;
    const double lagDenominator = lagTangent.squaredNorm();
    const Eigen::Matrix<double,6,6> lagProjection =
        lagTangent * lagTangent.transpose() / lagDenominator;
    const Eigen::Matrix<double,1,6> lagRow =
        lagTangent.transpose() / std::sqrt(lagDenominator);
    const Eigen::Matrix<double,6,6> scalarLagCost =
        10.0 * lagRow.transpose() * lagRow;
    const Eigen::Matrix<double,6,6> matrixLagCost =
        lagProjection.transpose()
        * (10.0 * Eigen::Matrix<double,6,6>::Identity())
        * lagProjection;
    if((scalarLagCost - matrixLagCost).norm() > 1e-12)
    {
        std::cerr << "Matrix lag cost does not reproduce the scalar isotropic cost.\n";
        return 7;
    }

    return 0;
}
