/** @file MpccQpConstraints.cpp
 * Unit tests for the fixed-variable box split used by the MPCC condensed QP.
 */

#include "MpccQpConstraints.h"

#include <cmath>
#include <iostream>
#include <string>

namespace {

using RobotLibrary::Control::detail::PinnedBoxSplit;
using RobotLibrary::Control::detail::split_pinned_box;

int failures = 0;

void check(const bool condition, const std::string &message)
{
    if(not condition)
    {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

Eigen::VectorXd vec(const std::initializer_list<double> values)
{
    Eigen::VectorXd result(static_cast<Eigen::Index>(values.size()));
    Eigen::Index i = 0;
    for(const double value : values) result(i++) = value;
    return result;
}

/** A pinned variable becomes one equality; free variables keep both bounds. */
void mixed_box_splits_by_variable()
{
    const Eigen::VectorXd lower = vec({-1.0, 0.5, -2.0});
    const Eigen::VectorXd upper = vec({ 1.0, 0.5,  3.0});
    const PinnedBoxSplit split = split_pinned_box(lower, upper);

    check(split.pinnedIndices.size() == 1 && split.pinnedIndices[0] == 1,
          "only the coincident-bound variable is pinned");
    check(split.equalityMatrix.rows() == 1 && split.equalityMatrix.cols() == 3,
          "one equality row spanning every variable");
    check(split.equalityMatrix(0, 1) == 1.0
          && split.equalityMatrix(0, 0) == 0.0 && split.equalityMatrix(0, 2) == 0.0,
          "the equality row is the pinned variable's unit vector");
    check(std::abs(split.equalityVector(0) - 0.5) <= 1e-15,
          "the equality fixes the variable to its coincident bound");
    check(split.inequalityMatrix.rows() == 4,
          "two free variables leave four inequality rows");

    // Every free row must still bound exactly what the original box did.
    const Eigen::VectorXd inside = vec({0.25, 0.5, 2.0});
    check((split.inequalityMatrix * inside - split.inequalityVector).maxCoeff() <= 0.0,
          "an interior point satisfies the retained inequalities");
    const Eigen::VectorXd outside = vec({1.5, 0.5, 2.0});
    check((split.inequalityMatrix * outside - split.inequalityVector).maxCoeff() > 0.0,
          "a point outside a free bound still violates a retained inequality");
}

/** With no coincident bounds the split must reproduce the plain box. */
void free_box_is_untouched()
{
    const Eigen::VectorXd lower = vec({-1.0, -2.0});
    const Eigen::VectorXd upper = vec({ 1.0,  2.0});
    const PinnedBoxSplit split = split_pinned_box(lower, upper);

    check(split.equalityMatrix.rows() == 0 && split.equalityVector.size() == 0,
          "no equalities when every variable has a genuine interval");
    check(split.inequalityMatrix.rows() == 4, "all four box rows are retained");
}

/**
 * The near-end-of-path case: the MPCC lower progress bound drops below the
 * upper one, and such a variable must NOT be pinned. Fixing a variable that
 * still has a real interval would trade a conditioning defect for a
 * correctness defect.
 */
void nearly_equal_bounds_stay_free_by_default()
{
    const Eigen::VectorXd lower = vec({1.0 - 1e-9});
    const Eigen::VectorXd upper = vec({1.0});
    const PinnedBoxSplit split = split_pinned_box(lower, upper);

    check(split.equalityMatrix.rows() == 0,
          "a narrow but nonzero interval is not pinned at the default tolerance");
    check(split.inequalityMatrix.rows() == 2, "both bounds are retained");

    const PinnedBoxSplit tolerant = split_pinned_box(lower, upper, 1e-6);
    check(tolerant.equalityMatrix.rows() == 1,
          "an explicit tolerance may pin a narrow interval");
    check(std::abs(tolerant.equalityVector(0) - 0.5 * (1.0 + 1.0 - 1e-9)) <= 1e-15,
          "a tolerated pin fixes the interval midpoint");
}

/** An inverted box is a configuration error, not something to silently pin. */
void inverted_bounds_are_rejected()
{
    bool threw = false;
    try { (void)split_pinned_box(vec({1.0}), vec({-1.0})); }
    catch(const std::invalid_argument &) { threw = true; }
    check(threw, "a lower bound above its upper bound is rejected");

    threw = false;
    try { (void)split_pinned_box(vec({0.0, 0.0}), vec({1.0})); }
    catch(const std::invalid_argument &) { threw = true; }
    check(threw, "mismatched bound sizes are rejected");

    threw = false;
    try
    {
        (void)split_pinned_box(vec({std::numeric_limits<double>::quiet_NaN()}), vec({1.0}));
    }
    catch(const std::invalid_argument &) { threw = true; }
    check(threw, "non-finite bounds are rejected");
}

/** The split describes the same feasible set as the box it replaces. */
void split_preserves_the_feasible_set()
{
    const Eigen::VectorXd lower = vec({-1.0, 2.0, 0.0});
    const Eigen::VectorXd upper = vec({ 1.0, 2.0, 4.0});
    const PinnedBoxSplit split = split_pinned_box(lower, upper);

    const auto feasibleInSplit = [&](const Eigen::VectorXd &x)
    {
        const bool equalitiesHold = split.equalityMatrix.rows() == 0
            || (split.equalityMatrix * x - split.equalityVector).cwiseAbs().maxCoeff() <= 1e-12;
        return equalitiesHold
            && (split.inequalityMatrix * x - split.inequalityVector).maxCoeff() <= 1e-12;
    };
    const auto feasibleInBox = [&](const Eigen::VectorXd &x)
    {
        return (x - lower).minCoeff() >= -1e-12 && (upper - x).minCoeff() >= -1e-12;
    };

    const Eigen::VectorXd points[] = {
        vec({ 0.0, 2.0, 2.0}),   // interior
        vec({ 1.0, 2.0, 4.0}),   // on the free upper bounds
        vec({-1.0, 2.0, 0.0}),   // on the free lower bounds
        vec({ 0.0, 2.5, 2.0}),   // violates the pin
        vec({ 2.0, 2.0, 2.0}),   // violates a free bound
    };
    for(const Eigen::VectorXd &x : points)
    {
        check(feasibleInSplit(x) == feasibleInBox(x),
              "split and box agree on feasibility");
    }
}

} // namespace

int main()
{
    mixed_box_splits_by_variable();
    free_box_is_untouched();
    nearly_equal_bounds_stay_free_by_default();
    inverted_bounds_are_rejected();
    split_preserves_the_feasible_set();

    if(failures == 0) std::cout << "mpcc_qp_constraints_test: all checks passed\n";
    return failures == 0 ? 0 : 1;
}
