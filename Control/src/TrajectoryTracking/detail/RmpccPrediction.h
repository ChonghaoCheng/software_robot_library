/**
 * @file TrajectoryTracking/detail/RmpccPrediction.h
 * @brief Exact discrete SE(3) rollout and local RTI linearisation for RMPCC.
 */

#ifndef RMPCC_PREDICTION_H
#define RMPCC_PREDICTION_H

#include <Math/MathFunctions.h>

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace RobotLibrary { namespace Control {

using RmpccStateVector = Eigen::Matrix<double,7,1>;  // [log(E); s]
using RmpccInputVector = Eigen::Matrix<double,7,1>;  // [u_body; sdot]

struct RmpccStageLinearization
{
    RmpccStateVector nominalNext = RmpccStateVector::Zero();
    Eigen::Matrix<double,7,7> stateJacobian = Eigen::Matrix<double,7,7>::Zero();
    Eigen::Matrix<double,7,7> inputJacobian = Eigen::Matrix<double,7,7>::Zero();
};

inline std::pair<double,double>
rmpcc_progress_rate_perturbations(const double progress,
                                  const double progressRate,
                                  const double dt,
                                  const double finiteDifferenceStep)
{
    const double boundedProgress = std::clamp(progress, 0.0, 1.0);
    const double minimumRate = -boundedProgress / dt;
    const double maximumRate = (1.0 - boundedProgress) / dt;
    return {
        std::clamp(progressRate - finiteDifferenceStep, minimumRate, maximumRate),
        std::clamp(progressRate + finiteDifferenceStep, minimumRate, maximumRate)
    };
}

/**
 * @brief Express the reference body tangent in the current endpoint body frame.
 *
 * For E = T_ref^-1 T and reference body tangent tau, keeping E constant
 * requires u = Ad(E^-1) tau sdot.
 */
inline Eigen::Matrix<double,6,1>
rmpcc_transport_reference_tangent(
    const Eigen::Matrix<double,6,1> &error,
    const Eigen::Matrix<double,6,1> &referenceTangent)
{
    const Eigen::Matrix4d relative =
        RobotLibrary::Math::se3_exponential(error);
    return RobotLibrary::Math::adjoint(
               RobotLibrary::Math::se3_inverse(relative))
           * referenceTangent;
}

/**
 * @brief Positive path tangent in logarithmic error coordinates.
 *
 * At fixed endpoint pose, d log(T_ref(s)^-1 T) / ds is the negative of this
 * vector. The positive form preserves the convention that positive lag means
 * the endpoint lies ahead of the reference:
 *
 *   g(e,s) = J_r(e)^-1 Ad(E^-1) tau(s).
 */
inline Eigen::Matrix<double,6,1>
rmpcc_error_coordinate_path_tangent(
    const Eigen::Matrix<double,6,1> &error,
    const Eigen::Matrix<double,6,1> &referenceTangent)
{
    // J_r(e)^-1 Ad(Exp(-e)) = J_l(e)^-1.  In the controller's local
    // operating region, evaluate the inverse left Jacobian directly with the
    // Bernoulli series.  This preserves the complete state dependence while
    // avoiding one SE(3) exponential and one 6x6 FullPivLU per residual call.
    // Retain the established general implementation outside that region.
    if(error.tail<3>().norm() <= 0.5)
    {
        using Matrix6d = Eigen::Matrix<double,6,6>;
        const Matrix6d ad = RobotLibrary::Math::se3_adjoint_matrix(error);
        const Matrix6d ad2 = ad * ad;
        const Matrix6d ad4 = ad2 * ad2;
        const Matrix6d ad6 = ad4 * ad2;
        const Matrix6d ad8 = ad4 * ad4;
        const Matrix6d ad10 = ad8 * ad2;
        const Matrix6d ad12 = ad10 * ad2;
        const Matrix6d inverseLeftJacobian =
            Matrix6d::Identity() - 0.5 * ad + (1.0 / 12.0) * ad2
            - (1.0 / 720.0) * ad4 + (1.0 / 30240.0) * ad6
            - (1.0 / 1209600.0) * ad8 + (1.0 / 47900160.0) * ad10
            - (691.0 / 1307674368000.0) * ad12;
        return inverseLeftJacobian * referenceTangent;
    }
    return RobotLibrary::Math::se3_right_jacobian_inverse(error)
           * rmpcc_transport_reference_tangent(error, referenceTangent);
}

/**
 * @brief Additive forward-Euler approximation of the same logarithmic error
 *        state used by the exact SE(3) rollout.
 *
 * For E = T_ref(s)^-1 T and e = Log(E), the local model is
 *
 *   e_next = e + dt (J_r(e)^-1 u - g(e,s) sdot),
 *   g(e,s) = J_r(e)^-1 Ad(E^-1) tau(s).
 *
 * This deliberately replaces group multiplication by vector addition while
 * preserving the state, input, tangent, and progress conventions.  It exists
 * for the predictor x cost-geometry factorial experiment; production RMPCC
 * continues to default to the exact rollout.
 */
template<typename ReferenceTangentFunction>
RmpccStateVector
rmpcc_additive_state_step(const RmpccStateVector &state,
                          const RmpccInputVector &input,
                          const double dt,
                          ReferenceTangentFunction &&referenceTangent)
{
    if(not state.allFinite() or not input.allFinite()
       or not std::isfinite(dt) or dt <= 0.0)
    {
        throw std::invalid_argument(
            "[ERROR] [RMPCC PREDICTION] State, input and positive dt must be finite.");
    }

    const double progress = std::clamp(state(6), 0.0, 1.0);
    const Eigen::Matrix<double,6,6> rightJacobianInverse =
        RobotLibrary::Math::se3_right_jacobian_inverse(state.head<6>());
    const Eigen::Matrix<double,6,1> tangent = referenceTangent(progress);

    RmpccStateVector nextState;
    nextState.head<6>() = state.head<6>() + dt * (
        rightJacobianInverse * input.head<6>()
        - rmpcc_error_coordinate_path_tangent(state.head<6>(), tangent) * input(6));
    nextState(6) = std::clamp(progress + dt * input(6), 0.0, 1.0);
    return nextState;
}

/** Central-difference RTI linearisation of rmpcc_additive_state_step(). */
template<typename ReferenceTangentFunction>
RmpccStageLinearization
rmpcc_linearize_additive_state_step(
    const RmpccStateVector &state,
    const RmpccInputVector &input,
    const double dt,
    const double finiteDifferenceStep,
    ReferenceTangentFunction &&referenceTangent)
{
    if(not std::isfinite(finiteDifferenceStep) or finiteDifferenceStep <= 0.0)
    {
        throw std::invalid_argument(
            "[ERROR] [RMPCC PREDICTION] finiteDifferenceStep must be positive.");
    }

    auto &&tangent = referenceTangent;
    RmpccStageLinearization result;
    result.nominalNext = rmpcc_additive_state_step(state, input, dt, tangent);
    for(int column = 0; column < 7; ++column)
    {
        RmpccStateVector plus = state;
        RmpccStateVector minus = state;
        plus(column) += finiteDifferenceStep;
        minus(column) -= finiteDifferenceStep;
        if(column == 6)
        {
            plus(column) = std::clamp(plus(column), 0.0, 1.0);
            minus(column) = std::clamp(minus(column), 0.0, 1.0);
        }
        const double denominator = plus(column) - minus(column);
        if(std::abs(denominator) > 1e-15)
        {
            result.stateJacobian.col(column) =
                (rmpcc_additive_state_step(plus, input, dt, tangent)
                 - rmpcc_additive_state_step(minus, input, dt, tangent)) / denominator;
        }
    }
    for(int column = 0; column < 7; ++column)
    {
        RmpccInputVector plus = input;
        RmpccInputVector minus = input;
        plus(column) += finiteDifferenceStep;
        minus(column) -= finiteDifferenceStep;
        double denominator = 2.0 * finiteDifferenceStep;
        if(column == 6)
        {
            const auto perturbations = rmpcc_progress_rate_perturbations(
                state(6), input(6), dt, finiteDifferenceStep);
            minus(6) = perturbations.first;
            plus(6) = perturbations.second;
            denominator = plus(6) - minus(6);
        }
        if(std::abs(denominator) > 1e-15)
        {
            result.inputJacobian.col(column) =
                (rmpcc_additive_state_step(state, plus, dt, tangent)
                 - rmpcc_additive_state_step(state, minus, dt, tangent))
                / denominator;
        }
    }
    return result;
}

/**
 * @brief Propagate one RMPCC stage by direct Lie-group multiplication.
 *
 * With E_k = T_ref(s_k)^-1 T_k and a piecewise-constant current-endpoint
 * body twist u_k,
 *
 *   E_(k+1) = T_ref(s_(k+1))^-1 T_ref(s_k)
 *             E_k Exp(dt u_k),
 *   s_(k+1) = clamp(s_k + dt sdot_k, 0, 1).
 *
 * The returned state stores log(E_(k+1)); the nominal pose rollout itself is
 * performed on SE(3), not by integrating logarithm coordinates.
 */
template<typename ReferenceTransformFunction>
RmpccStateVector
rmpcc_exact_state_step(const RmpccStateVector &state,
                       const RmpccInputVector &input,
                       const double dt,
                       ReferenceTransformFunction &&referenceTransform)
{
    if(not state.allFinite() or not input.allFinite()
       or not std::isfinite(dt) or dt <= 0.0)
    {
        throw std::invalid_argument(
            "[ERROR] [RMPCC PREDICTION] State, input and positive dt must be finite.");
    }

    const double progress = std::clamp(state(6), 0.0, 1.0);
    const double nextProgress = std::clamp(progress + dt * input(6), 0.0, 1.0);
    const Eigen::Matrix4d currentReference = referenceTransform(progress);
    const Eigen::Matrix4d nextReference = referenceTransform(nextProgress);
    const Eigen::Matrix4d relativeTransform =
        RobotLibrary::Math::se3_exponential(state.head<6>());
    const Eigen::Matrix4d nextRelativeTransform =
        RobotLibrary::Math::se3_inverse(nextReference)
        * currentReference
        * relativeTransform
        * RobotLibrary::Math::se3_exponential(dt * input.head<6>());

    RmpccStateVector nextState;
    nextState.head<6>() = RobotLibrary::Math::se3_logarithm(nextRelativeTransform);
    nextState(6) = nextProgress;
    return nextState;
}

/**
 * @brief Central-difference linearisation of the exact SE(3) stage map.
 *
 * This linearisation is used only for the RTI correction QP. The nominal
 * horizon is always generated by rmpcc_exact_state_step().
 */
template<typename ReferenceTransformFunction>
RmpccStageLinearization
rmpcc_linearize_exact_state_step(const RmpccStateVector &state,
                                 const RmpccInputVector &input,
                                 const double dt,
                                 const double finiteDifferenceStep,
                                 ReferenceTransformFunction &&referenceTransform)
{
    if(not std::isfinite(finiteDifferenceStep) or finiteDifferenceStep <= 0.0)
    {
        throw std::invalid_argument(
            "[ERROR] [RMPCC PREDICTION] finiteDifferenceStep must be positive.");
    }

    auto &&reference = referenceTransform;
    RmpccStageLinearization result;
    result.nominalNext = rmpcc_exact_state_step(state, input, dt, reference);

    for(int column = 0; column < 7; ++column)
    {
        RmpccStateVector plus = state;
        RmpccStateVector minus = state;
        plus(column) += finiteDifferenceStep;
        minus(column) -= finiteDifferenceStep;
        if(column == 6)
        {
            plus(column) = std::clamp(plus(column), 0.0, 1.0);
            minus(column) = std::clamp(minus(column), 0.0, 1.0);
        }
        const double denominator = plus(column) - minus(column);
        if(std::abs(denominator) > 1e-15)
        {
            result.stateJacobian.col(column) =
                (rmpcc_exact_state_step(plus, input, dt, reference)
                 - rmpcc_exact_state_step(minus, input, dt, reference)) / denominator;
        }
    }

    for(int column = 0; column < 7; ++column)
    {
        RmpccInputVector plus = input;
        RmpccInputVector minus = input;
        plus(column) += finiteDifferenceStep;
        minus(column) -= finiteDifferenceStep;
        double denominator = 2.0 * finiteDifferenceStep;
        if(column == 6)
        {
            const auto perturbations = rmpcc_progress_rate_perturbations(
                state(6), input(6), dt, finiteDifferenceStep);
            minus(6) = perturbations.first;
            plus(6) = perturbations.second;
            denominator = plus(6) - minus(6);
        }
        if(std::abs(denominator) > 1e-15)
        {
            result.inputJacobian.col(column) =
                (rmpcc_exact_state_step(state, plus, dt, reference)
                 - rmpcc_exact_state_step(state, minus, dt, reference))
                / denominator;
        }
    }

    return result;
}

} } // namespace RobotLibrary::Control

#endif
