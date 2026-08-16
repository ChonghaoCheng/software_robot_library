#include "RmpccQpConstraints.h"

#include <Math/QPSolver.h>

#include <Eigen/Core>
#include <Eigen/LU>

#include <cmath>
#include <iostream>
#include <stdexcept>

int main()
{
    using RobotLibrary::Control::RmpccQpConstraints;
    using RobotLibrary::Control::rmpcc_build_qp_constraints;
    using RobotLibrary::Control::rmpcc_fixed_progress_schedule_satisfies_limits;

    constexpr int horizon = 3;
    constexpr int dim = 7 * horizon;
    constexpr int progressOffset = 6 * horizon;
    constexpr double dt = 0.002;
    Eigen::VectorXd lower = Eigen::VectorXd::Constant(dim, -0.2);
    Eigen::VectorXd upper = Eigen::VectorXd::Constant(dim, 0.2);
    const Eigen::Vector3d fixedRates(0.05, 0.05, 0.05);
    lower.tail(horizon) = fixedRates;
    upper.tail(horizon) = fixedRates;

    const RmpccQpConstraints fixed = rmpcc_build_qp_constraints(
        lower, upper, horizon, dt, 0.9, 0.1, true, fixedRates);
    if(fixed.Aeq.rows() != horizon
       or Eigen::FullPivLU<Eigen::MatrixXd>(fixed.Aeq).rank() != horizon)
    {
        std::cerr << "Fixed-progress equalities are not independent.\n";
        return 1;
    }
    Eigen::VectorXd nominal = Eigen::VectorXd::Zero(dim);
    nominal.tail(horizon) = fixedRates;
    if((fixed.Aeq * nominal - fixed.yeq).cwiseAbs().maxCoeff() > 1e-14)
    {
        std::cerr << "Nominal point does not satisfy fixed-progress equalities.\n";
        return 2;
    }
    if(fixed.Bineq.rows() != 2 * progressOffset
       or fixed.Bineq.rightCols(horizon).cwiseAbs().maxCoeff() != 0.0)
    {
        std::cerr << "Fixed-progress bound inequalities were not removed.\n";
        return 3;
    }
    if(not rmpcc_fixed_progress_schedule_satisfies_limits(
           fixedRates, dt, 0.9, 0.1)
       or rmpcc_fixed_progress_schedule_satisfies_limits(
           fixedRates, dt, 1e-5, 0.1)
       or rmpcc_fixed_progress_schedule_satisfies_limits(
           fixedRates, dt, 0.9, 1e-5))
    {
        std::cerr << "Fixed-progress implied-limit validation is incorrect.\n";
        return 4;
    }

    const RmpccQpConstraints optimized = rmpcc_build_qp_constraints(
        lower, upper, horizon, dt, 0.9, 0.1, false, fixedRates);
    Eigen::MatrixXd expectedB = Eigen::MatrixXd::Zero(2 * dim + 2, dim);
    Eigen::VectorXd expectedZ = Eigen::VectorXd::Zero(2 * dim + 2);
    expectedB.topRows(dim).setIdentity();
    expectedZ.head(dim) = upper;
    expectedB.block(dim, 0, dim, dim) = -Eigen::MatrixXd::Identity(dim, dim);
    expectedZ.segment(dim, dim) = -lower;
    expectedB.block(2 * dim, progressOffset, 1, horizon).setConstant(dt);
    expectedB(2 * dim + 1, progressOffset) = dt;
    expectedZ(2 * dim) = 0.9;
    expectedZ(2 * dim + 1) = 0.1;
    if(optimized.Aeq.rows() != 0 or optimized.yeq.size() != 0
       or (optimized.Bineq - expectedB).norm() != 0.0
       or (optimized.zineq - expectedZ).norm() != 0.0)
    {
        std::cerr << "Optimized-progress inequality representation changed.\n";
        return 5;
    }

    // Reduced rotation-rich-like working set: linear and angular twist start
    // at their bounds while progress is fixed. Encoding fixed progress as two
    // inequalities creates four active rows for only three variables.
    Eigen::Matrix3d H;
    H << 2.0, 0.2, 0.1,
         0.2, 1.5, -0.1,
         0.1, -0.1, 1.0;
    Eigen::Vector3d f(-1.0, 1.0, -0.02);
    Eigen::Vector3d x0(0.2, -0.2, 0.05);
    Eigen::Matrix<double,6,3> oldB;
    oldB.topRows<3>().setIdentity();
    oldB.bottomRows<3>() = -Eigen::Matrix3d::Identity();
    Eigen::Vector<double,6> oldZ;
    oldZ << 0.2, 0.2, 0.05, 0.2, 0.2, -0.05;
    SolverOptions<double> options;
    options.method = "active set";
    options.stepSizeTolerance = 1e-5;
    options.maxSteps = 50;
    QPSolver<double> solver(options);
    bool oldRepresentationFailed = false;
    try
    {
        (void)solver.solve(H, f, oldB, oldZ, x0);
    }
    catch(const std::invalid_argument &)
    {
        oldRepresentationFailed = true;
    }
    if(not oldRepresentationFailed)
    {
        std::cerr << "Reduced paired-inequality regression no longer reproduces failure.\n";
        return 6;
    }

    Eigen::Matrix<double,1,3> Aeq;
    Aeq << 0.0, 0.0, 1.0;
    Eigen::Vector<double,1> yeq;
    yeq << 0.05;
    Eigen::Matrix<double,4,3> twistB;
    twistB.topRows<2>() = oldB.topRows<2>();
    twistB.bottomRows<2>() = oldB.middleRows<2>(3);
    Eigen::Vector4d twistZ(0.2, 0.2, 0.2, 0.2);
    const Eigen::Vector3d corrected =
        solver.solve(H, f, Aeq, yeq, twistB, twistZ, x0);
    const double equalityResidual = std::abs((Aeq * corrected - yeq)(0));
    const double inequalityViolation = (twistB * corrected - twistZ).maxCoeff();
    const double objective = 0.5 * corrected.dot(H * corrected) + f.dot(corrected);
    if(not corrected.allFinite() or not std::isfinite(objective)
       or equalityResidual >= 1e-9 or inequalityViolation >= 1e-9)
    {
        std::cerr << "Equality-form reduced regression is invalid: eq="
                  << equalityResidual << " ineq=" << inequalityViolation << '\n';
        return 7;
    }
    return 0;
}
