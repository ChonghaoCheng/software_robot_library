#include "RmpccReferenceMotion.h"

#include <Math/MathFunctions.h>

#include <Eigen/Core>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

namespace {

using Matrix4 = Eigen::Matrix4d;
using Vector6 = Eigen::Vector<double,6>;
using RobotLibrary::Control::RmpccInputVector;
using RobotLibrary::Control::RmpccStateVector;

Matrix4 exp(const Vector6 &value)
{
    return RobotLibrary::Math::se3_exponential(value);
}

Matrix4 inverse(const Matrix4 &value)
{
    return RobotLibrary::Math::se3_inverse(value);
}

Vector6 log(const Matrix4 &value)
{
    return RobotLibrary::Math::se3_logarithm(value);
}

double metric_norm(const Vector6 &value)
{
    return std::sqrt(value.head<3>().squaredNorm()
                     + 0.04 * value.tail<3>().squaredNorm());
}

} // namespace

int main()
{
    using namespace RobotLibrary::Control;
    Vector6 pathLinear;
    pathLinear << 0.18, -0.07, 0.11, 0.0, 0.0, 0.0;
    Vector6 pathRich;
    pathRich << -0.03, 0.09, 0.04, 0.62, -0.41, 0.33;
    Vector6 pathCurve;
    pathCurve << 0.02, -0.04, 0.03, -0.31, 0.48, 0.27;

    const auto makeReference = [&](const bool rich)
    {
        return [=](const double progress) -> Matrix4
        {
            const double s = std::clamp(progress, 0.0, 1.0);
            return (exp(pathLinear * s)
                   * exp((rich ? pathRich : Vector6::Zero()) * s)
                   * exp((rich ? pathCurve : Vector6::Zero()) * (s * s))).eval();
        };
    };

    const std::array<Vector6,4> parentSteps{{
        Vector6::Zero(),
        (Vector6() << 0.001, -0.0004, 0.0002, 0, 0, 0).finished(),
        (Vector6() << 0, 0, 0, 0.002, -0.001, 0.0015).finished(),
        (Vector6() << -0.0005, 0.0008, 0.0003, -0.001, 0.0017, 0.0006).finished()
    }};
    Vector6 parentOrigin;
    parentOrigin << 0.2, -0.1, 0.3, 0.15, -0.08, 0.11;

    constexpr double dt = 0.002;
    constexpr double progress = 0.41;
    constexpr double rate = 0.027;
    double maxFactorization = 0.0;
    double maxStaticParity = 0.0;
    double maxInvariant = 0.0;
    double maxLinearizationRelative = 0.0;

    const std::array<Vector6,4> errors{{
        Vector6::Zero(),
        (Vector6() << 0.001, -0.0004, 0.0007, 0, 0, 0).finished(),
        (Vector6() << 0, 0, 0, 0.005, -0.003, 0.004).finished(),
        (Vector6() << 0.0006, -0.0008, 0.0003, 0.004, 0.002, -0.003).finished()
    }};

    for(const bool rich : {false, true})
    {
        const auto reference = makeReference(rich);
        for(const Vector6 &parentStep : parentSteps)
        {
            const auto parent = [&](const int stage) -> Matrix4
            {
                return (exp(parentOrigin)
                        * exp(parentStep * static_cast<double>(stage))).eval();
            };
            const double next = std::clamp(progress + dt * rate, 0.0, 1.0);
            const Matrix4 direct = inverse(reference(progress))
                                   * inverse(parent(0)) * parent(1) * reference(next);
            const Matrix4 candidate = rmpcc_active_reference_displacement(
                progress, rate, dt, 0, reference, parent);
            maxFactorization = std::max(
                maxFactorization, log(inverse(direct) * candidate).norm());

            for(const Vector6 &error : errors)
            {
                RmpccStateVector state = RmpccStateVector::Zero();
                state.head<6>() = error;
                state(6) = progress;
                RmpccInputVector input = RmpccInputVector::Zero();
                input.head<6>() = rmpcc_finite_stage_feedforward(
                    error, progress, rate, dt, 0, reference, parent);
                input(6) = rate;
                const RmpccStateVector nextState = rmpcc_finite_stage_state_step(
                    state, input, dt, 0, reference, parent);
                maxInvariant = std::max(
                    maxInvariant,
                    metric_norm(log(inverse(exp(error)) * exp(nextState.head<6>()))));
            }

            RmpccStateVector state;
            state << 0.012, -0.008, 0.006, 0.031, -0.024, 0.017, progress;
            RmpccInputVector input;
            input << 0.018, -0.011, 0.009, 0.041, -0.027, 0.019, rate;
            const auto linearization =
                rmpcc_linearize_finite_stage_path_velocity_residual(
                    state, input, dt, 0, 1e-6, reference, parent);
            RmpccStateVector stateDirection;
            stateDirection << 0.21, -0.37, 0.18, 0.29, -0.16, 0.31, 0.23;
            stateDirection.normalize();
            RmpccInputVector inputDirection;
            inputDirection << -0.14, 0.26, 0.33, -0.28, 0.19, 0.22, -0.41;
            inputDirection.normalize();
            const Vector6 predicted = linearization.stateJacobian * stateDirection
                                      + linearization.inputJacobian * inputDirection;
            double bestRelative = 1e100;
            for(const double outerStep : {2e-7, 5e-7, 2e-6})
            {
                const Vector6 measured =
                    (rmpcc_finite_stage_path_velocity_residual(
                         state + outerStep * stateDirection,
                         input + outerStep * inputDirection,
                         dt, 0, reference, parent)
                     - rmpcc_finite_stage_path_velocity_residual(
                         state - outerStep * stateDirection,
                         input - outerStep * inputDirection,
                         dt, 0, reference, parent)) / (2.0 * outerStep);
                bestRelative = std::min(
                    bestRelative,
                    (measured - predicted).norm() / std::max(1.0, measured.norm()));
            }
            maxLinearizationRelative = std::max(
                maxLinearizationRelative, bestRelative);
        }

        const auto identityParent = [](const int) { return Matrix4::Identity(); };
        RmpccStateVector state;
        state << 0.008, -0.005, 0.003, 0.012, -0.009, 0.007, progress;
        RmpccInputVector input;
        input << 0.015, -0.012, 0.01, 0.021, -0.016, 0.013, rate;
        const RmpccStateVector candidate = rmpcc_finite_stage_state_step(
            state, input, dt, 0, reference, identityParent);
        const RmpccStateVector existing = rmpcc_exact_state_step(
            state, input, dt, reference);
        maxStaticParity = std::max(maxStaticParity, (candidate - existing).norm());
    }

    std::cout.precision(17);
    std::cout << "factorization_max=" << maxFactorization << '\n'
              << "static_predictor_parity_max=" << maxStaticParity << '\n'
              << "constant_E_metric_max=" << maxInvariant << '\n'
              << "linearization_relative_max=" << maxLinearizationRelative << '\n';
    if(maxFactorization >= 1e-13 || maxStaticParity >= 1e-12
       || maxInvariant > 1e-12 || maxLinearizationRelative > 5e-5)
    {
        return 1;
    }
    return 0;
}
