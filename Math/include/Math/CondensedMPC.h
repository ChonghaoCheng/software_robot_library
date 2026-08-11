/**
 * @file    CondensedMPC.h
 * @author  (added by AI)
 * @brief   Reusable assembly helpers for condensed (dense) linear MPC problems.
 *
 * @details These free functions build the matrices that every condensed linear-MPC
 *          controller needs: the stacked prediction matrices for a linear model
 *          x_{k+1} = A x_k + B u_k, block-diagonal weight matrices, and box
 *          constraints in the B*z <= c form expected by QPSolver. They are pure
 *          Eigen and have no dependency on Model/Control, so any controller can
 *          share them instead of re-deriving the same algebra.
 *
 * @copyright (c) 2026
 *
 * @license   OSCL - Free for non-commercial open-source use only.
 */

#ifndef CONDENSED_MPC_H
#define CONDENSED_MPC_H

#include <Eigen/Dense>

namespace RobotLibrary { namespace Math {

/**
 * @brief Stacked prediction matrices for a condensed linear MPC.
 *
 * For x_{k+1} = A x_k + B u_k, the stacked future state
 * X = [x_1; x_2; ...; x_N] and stacked input U = [u_0; u_1; ...; u_{N-1}] satisfy
 *     X = stateTransition * x_0 + inputResponse * U.
 */
struct CondensedPrediction
{
    Eigen::MatrixXd stateTransition;   ///< Ax: (N*nx) x nx, with block k equal to A^{k+1}.
    Eigen::MatrixXd inputResponse;     ///< Bu: (N*nx) x (N*nu), block (k,j) = A^{k-j} B for j <= k, else 0.
};

/**
 * @brief Build the condensed prediction matrices for x_{k+1} = A x_k + B u_k.
 * @param A The nx x nx state transition matrix.
 * @param B The nx x nu input matrix.
 * @param horizon The prediction horizon N (>= 1; values of 0 are treated as 1).
 * @return A CondensedPrediction holding stateTransition (Ax) and inputResponse (Bu).
 */
CondensedPrediction
condense_prediction(const Eigen::MatrixXd &A,
                    const Eigen::MatrixXd &B,
                    const unsigned int    &horizon);

/**
 * @brief Replicate a square block along the diagonal `count` times.
 * @param block The square block to repeat (e.g. a per-stage Q or R weight).
 * @param count How many times to place it on the diagonal.
 * @return A (count*n) x (count*n) block-diagonal matrix (zero off the blocks).
 */
Eigen::MatrixXd
block_diagonal(const Eigen::MatrixXd &block,
               const unsigned int    &count);

/**
 * @brief Box constraints lower <= z <= upper expressed as B*z <= c.
 *
 * With B = [I; -I] and c = [upper; -lower], the single inequality B*z <= c is
 * exactly the element-wise box, in the form expected by QPSolver::solve().
 */
struct BoxConstraint
{
    Eigen::MatrixXd constraintMatrix;  ///< B: (2n) x n, [I; -I].
    Eigen::VectorXd constraintVector;  ///< c: 2n, [upper; -lower].
};

/**
 * @brief Build the (B, c) pair representing lower <= z <= upper.
 * @param lower Element-wise lower bounds (length n).
 * @param upper Element-wise upper bounds (length n).
 * @return A BoxConstraint with constraintMatrix [I; -I] and constraintVector [upper; -lower].
 */
BoxConstraint
box_constraint(const Eigen::VectorXd &lower,
               const Eigen::VectorXd &upper);

} } // namespace

#endif
