/**
 * @file TrajectoryTracking/detail/LieAlgebraMPCPrediction.h
 * @brief Linearized SE(3) error dynamics used by Lie-algebra MPC.
 */

#ifndef LIE_ALGEBRA_MPC_PREDICTION_H
#define LIE_ALGEBRA_MPC_PREDICTION_H

#include <Math/MathFunctions.h>

#include <Eigen/Core>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace RobotLibrary { namespace Control {

using LieAlgebraMPCVector = Eigen::Matrix<double,6,1>;
using LieAlgebraMPCMatrix = Eigen::Matrix<double,6,6>;

/**
 * @brief Forward-Euler linearization of a left-invariant tracking error.
 *
 * For E = T_ref^-1 T, e = Log(E), body twists xi_ref and
 * xi = xi_ref + delta_xi, the first-order continuous dynamics are
 *
 *   e_dot = -ad(xi_ref) e + delta_xi.
 *
 * Thus e_next = A e + B delta_xi.
 */
inline std::pair<LieAlgebraMPCMatrix, LieAlgebraMPCMatrix>
lie_algebra_mpc_stage(const LieAlgebraMPCVector &referenceBodyTwist,
                      const double dt)
{
    if(not referenceBodyTwist.allFinite() or not std::isfinite(dt) or dt <= 0.0)
    {
        throw std::invalid_argument(
            "[ERROR] [LIE ALGEBRA MPC] Reference twist and positive dt must be finite.");
    }

    const LieAlgebraMPCMatrix A = LieAlgebraMPCMatrix::Identity()
        - dt * RobotLibrary::Math::se3_adjoint_matrix(referenceBodyTwist);
    const LieAlgebraMPCMatrix B = dt * LieAlgebraMPCMatrix::Identity();
    return {A, B};
}

/** Convert an endpoint-origin twist from base axes to body axes. */
inline LieAlgebraMPCVector
endpoint_twist_in_body(const Eigen::Matrix3d &bodyRotationInBase,
                       const LieAlgebraMPCVector &twistInBase)
{
    LieAlgebraMPCVector result;
    result.head<3>() = bodyRotationInBase.transpose() * twistInBase.head<3>();
    result.tail<3>() = bodyRotationInBase.transpose() * twistInBase.tail<3>();
    return result;
}

/** Convert an endpoint-origin twist from body axes to base axes. */
inline LieAlgebraMPCVector
endpoint_twist_in_base(const Eigen::Matrix3d &bodyRotationInBase,
                       const LieAlgebraMPCVector &twistInBody)
{
    LieAlgebraMPCVector result;
    result.head<3>() = bodyRotationInBase * twistInBody.head<3>();
    result.tail<3>() = bodyRotationInBase * twistInBody.tail<3>();
    return result;
}

} } // namespace RobotLibrary::Control

#endif
