/**
 * @file MpccQpConstraints.h
 * @brief Box-constraint assembly helpers shared by the MPCC condensed QP.
 */

#ifndef CONTROL_TRAJECTORY_TRACKING_DETAIL_MPCC_QP_CONSTRAINTS_H
#define CONTROL_TRAJECTORY_TRACKING_DETAIL_MPCC_QP_CONSTRAINTS_H

#include <Eigen/Core>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace RobotLibrary { namespace Control { namespace detail {

/**
 * @brief A box `lower <= x <= upper` split into equality and inequality parts.
 *
 * Variables whose bounds coincide are returned as one independent equality row
 * each; every other variable keeps the usual opposing pair of inequality rows.
 */
struct PinnedBoxSplit
{
    Eigen::MatrixXd equalityMatrix;     ///< p x n, one unit row per pinned variable
    Eigen::VectorXd equalityVector;     ///< p, the value each pinned variable is fixed to
    Eigen::MatrixXd inequalityMatrix;   ///< 2*(n-p) x n, stacked as [I_free; -I_free]
    Eigen::VectorXd inequalityVector;   ///< 2*(n-p)
    std::vector<int> pinnedIndices;     ///< Ascending variable indices carried by equalityMatrix
};

/**
 * @brief Split a box constraint so that fixed variables become equalities.
 *
 * @details An active-set QP that represents a fixed variable as two opposing
 *          inequalities starts on a degenerate vertex: both rows are tight at
 *          any feasible seed, so the ratio test returns a zero step and each
 *          iteration is spent identifying one redundant row instead of making
 *          primal progress. Emitting one independent equality per fixed
 *          variable removes the mutually opposite active rows without changing
 *          the feasible set.
 *
 * @param lower     Lower bounds (n).
 * @param upper     Upper bounds (n).
 * @param tolerance Bounds closer than this are treated as pinned. The default
 *                  pins only exactly coincident bounds.
 * @return The equality/inequality split.
 * @throws std::invalid_argument on size mismatch, non-finite bounds, or a
 *         lower bound above its upper bound.
 */
inline PinnedBoxSplit
split_pinned_box(const Eigen::VectorXd &lower,
                 const Eigen::VectorXd &upper,
                 const double tolerance = 0.0)
{
    if(lower.size() != upper.size())
    {
        throw std::invalid_argument(
            "[ERROR] [MPCC QP CONSTRAINTS] split_pinned_box(): "
            "Lower and upper bound vectors must have the same size.");
    }
    if(not lower.allFinite() or not upper.allFinite())
    {
        throw std::invalid_argument(
            "[ERROR] [MPCC QP CONSTRAINTS] split_pinned_box(): "
            "Bounds must be finite.");
    }
    if(not std::isfinite(tolerance) or tolerance < 0.0)
    {
        throw std::invalid_argument(
            "[ERROR] [MPCC QP CONSTRAINTS] split_pinned_box(): "
            "Tolerance must be finite and nonnegative.");
    }

    const int n = static_cast<int>(lower.size());

    PinnedBoxSplit split;
    std::vector<int> freeIndices;
    freeIndices.reserve(static_cast<std::size_t>(n));

    for(int i = 0; i < n; ++i)
    {
        const double width = upper(i) - lower(i);
        if(width < -tolerance)
        {
            throw std::invalid_argument(
                "[ERROR] [MPCC QP CONSTRAINTS] split_pinned_box(): "
                "Lower bound exceeds upper bound for variable "
                + std::to_string(i) + " (" + std::to_string(lower(i))
                + " > " + std::to_string(upper(i)) + ").");
        }
        if(width <= tolerance) split.pinnedIndices.push_back(i);
        else                   freeIndices.push_back(i);
    }

    const int pinned = static_cast<int>(split.pinnedIndices.size());
    const int free = static_cast<int>(freeIndices.size());

    split.equalityMatrix = Eigen::MatrixXd::Zero(pinned, n);
    split.equalityVector = Eigen::VectorXd::Zero(pinned);
    for(int k = 0; k < pinned; ++k)
    {
        const int index = split.pinnedIndices[static_cast<std::size_t>(k)];
        split.equalityMatrix(k, index) = 1.0;
        // Exact for coincident bounds, and the interior point otherwise.
        split.equalityVector(k) = 0.5 * (lower(index) + upper(index));
    }

    split.inequalityMatrix = Eigen::MatrixXd::Zero(2 * free, n);
    split.inequalityVector = Eigen::VectorXd::Zero(2 * free);
    for(int k = 0; k < free; ++k)
    {
        const int index = freeIndices[static_cast<std::size_t>(k)];
        split.inequalityMatrix(k, index) = 1.0;                 //  x_i <=  upper_i
        split.inequalityVector(k) = upper(index);
        split.inequalityMatrix(free + k, index) = -1.0;          // -x_i <= -lower_i
        split.inequalityVector(free + k) = -lower(index);
    }

    return split;
}

} } } // namespace RobotLibrary::Control::detail

#endif
