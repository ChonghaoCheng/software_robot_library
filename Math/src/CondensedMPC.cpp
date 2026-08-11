/**
 * @file    CondensedMPC.cpp
 * @author  (added by AI)
 * @brief   Implementation of the condensed linear MPC assembly helpers.
 *
 * @copyright (c) 2026
 *
 * @license   OSCL - Free for non-commercial open-source use only.
 */

#include <Math/CondensedMPC.h>

#include <stdexcept>
#include <vector>

namespace RobotLibrary { namespace Math {

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                          Condensed prediction matrices  X = Ax*x0 + Bu*U                        //
////////////////////////////////////////////////////////////////////////////////////////////////////
CondensedPrediction
condense_prediction(const Eigen::MatrixXd &A,
                    const Eigen::MatrixXd &B,
                    const unsigned int    &horizon)
{
    const int nx = static_cast<int>(A.rows());
    const int nu = static_cast<int>(B.cols());

    if(A.rows() != A.cols())
    {
        throw std::invalid_argument("[ERROR] [CONDENSE PREDICTION] State matrix A must be square.");
    }
    if(B.rows() != A.rows())
    {
        throw std::invalid_argument("[ERROR] [CONDENSE PREDICTION] "
                                    "B must have the same number of rows as A.");
    }

    const int N = (horizon == 0) ? 1 : static_cast<int>(horizon);

    CondensedPrediction result;
    result.stateTransition = Eigen::MatrixXd::Zero(N * nx, nx);
    result.inputResponse   = Eigen::MatrixXd::Zero(N * nx, N * nu);

    // Precompute the powers A^0 .. A^N once, so each block is a lookup rather than a re-multiply.
    std::vector<Eigen::MatrixXd> power(static_cast<size_t>(N) + 1,
                                       Eigen::MatrixXd::Identity(nx, nx));
    for(int i = 1; i <= N; ++i)
    {
        power[static_cast<size_t>(i)] = A * power[static_cast<size_t>(i - 1)];
    }

    for(int k = 0; k < N; ++k)
    {
        result.stateTransition.block(k * nx, 0, nx, nx) = power[static_cast<size_t>(k + 1)];         // A^{k+1}

        for(int j = 0; j <= k; ++j)
        {
            result.inputResponse.block(k * nx, j * nu, nx, nu) =
                power[static_cast<size_t>(k - j)] * B;                                               // A^{k-j} B
        }
    }

    return result;
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                          Replicate a square block along the diagonal                           //
////////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::MatrixXd
block_diagonal(const Eigen::MatrixXd &block,
               const unsigned int    &count)
{
    if(block.rows() != block.cols())
    {
        throw std::invalid_argument("[ERROR] [BLOCK DIAGONAL] The block must be square.");
    }

    const int n = static_cast<int>(block.rows());
    const int N = (count == 0) ? 1 : static_cast<int>(count);

    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(N * n, N * n);
    for(int i = 0; i < N; ++i)
    {
        result.block(i * n, i * n, n, n) = block;
    }

    return result;
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                          Box constraints  lower <= z <= upper  as  B z <= c                     //
////////////////////////////////////////////////////////////////////////////////////////////////////
BoxConstraint
box_constraint(const Eigen::VectorXd &lower,
               const Eigen::VectorXd &upper)
{
    if(lower.size() != upper.size())
    {
        throw std::invalid_argument("[ERROR] [BOX CONSTRAINT] "
                                    "Lower and upper bound vectors must have the same size.");
    }

    const int n = static_cast<int>(lower.size());

    BoxConstraint result;
    result.constraintMatrix = Eigen::MatrixXd::Zero(2 * n, n);
    result.constraintMatrix.topRows(n).setIdentity();
    result.constraintMatrix.bottomRows(n) = -Eigen::MatrixXd::Identity(n, n);

    result.constraintVector.resize(2 * n);
    result.constraintVector.head(n) =  upper;
    result.constraintVector.tail(n) = -lower;

    return result;
}

} } // namespace
