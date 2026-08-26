#include <Math/QPSolver.h>

#include <Eigen/Core>

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

namespace {

bool near(const Eigen::VectorXd &a, const Eigen::VectorXd &b, double tol)
{
    return (a - b).norm() <= tol;
}

bool check_kkt(const Eigen::MatrixXd &H,
               const Eigen::VectorXd &f,
               const Eigen::MatrixXd &B,
               const Eigen::VectorXd &z,
               const Eigen::VectorXd &x,
               double tol)
{
    const Eigen::VectorXd slack = z - B * x;
    if(slack.minCoeff() < -tol) return false;

    std::vector<int> active;
    for(int i = 0; i < slack.size(); ++i)
    {
        if(slack(i) <= 5.0 * tol) active.push_back(i);
    }

    Eigen::VectorXd mu = Eigen::VectorXd::Zero(active.size());
    if(!active.empty())
    {
        Eigen::MatrixXd Ba(active.size(), x.size());
        for(int i = 0; i < static_cast<int>(active.size()); ++i)
        {
            Ba.row(i) = B.row(active[static_cast<size_t>(i)]);
        }
        mu = Ba.transpose().completeOrthogonalDecomposition().solve(-(H * x + f));
        if(mu.minCoeff() < -2e-6) return false;
        if((H * x + f + Ba.transpose() * mu).norm() > 1e-6) return false;
    }
    else if((H * x + f).norm() > 1e-6)
    {
        return false;
    }

    for(int i = 0; i < static_cast<int>(active.size()); ++i)
    {
        if(std::abs(mu(i) * slack(active[static_cast<size_t>(i)])) > 1e-7) return false;
    }
    return true;
}

QPSolver<double> solver()
{
    SolverOptions<double> options;
    options.method = "active set";
    options.stepSizeTolerance = 1e-10;
    options.maxSteps = 100;
    return QPSolver<double>(options);
}

bool unconstrained_pd_qp()
{
    Eigen::Matrix2d H;
    H << 4.0, 1.0,
         1.0, 3.0;
    Eigen::Vector2d f(-1.0, 2.0);
    const Eigen::Vector2d expected = -H.ldlt().solve(f);
    const Eigen::Vector2d actual = QPSolver<double>::solve(H, f);
    return near(actual, expected, 1e-10);
}

bool active_bound_optimum()
{
    Eigen::Matrix2d H = Eigen::Matrix2d::Identity();
    Eigen::Vector2d f(-2.0, 0.5);
    Eigen::Matrix<double, 1, 2> B;
    B << 1.0, 0.0;
    Eigen::Vector<double, 1> z;
    z << 0.25;
    Eigen::Vector2d x0(0.0, 0.0);
    const Eigen::Vector2d x = solver().solve(H, f, B, z, x0);
    Eigen::Vector2d expected(0.25, -0.5);
    return near(x, expected, 1e-8) && check_kkt(H, f, B, z, x, 1e-8);
}

bool releases_initially_active_constraint()
{
    Eigen::Matrix<double, 1, 1> H;
    H << 1.0;
    Eigen::Vector<double, 1> f;
    f << 0.0;
    Eigen::Matrix<double, 1, 1> B;
    B << 1.0;
    Eigen::Vector<double, 1> z;
    z << 0.5;
    Eigen::Vector<double, 1> x0;
    x0 << 0.5;
    const Eigen::Vector<double, 1> x = solver().solve(H, f, B, z, x0);
    return std::abs(x(0)) <= 1e-8 && check_kkt(H, f, B, z, x, 1e-8);
}

bool redundant_boundary_seed_is_well_posed()
{
    Eigen::Matrix2d H = Eigen::Matrix2d::Identity();
    Eigen::Vector2d f(-1.0, 0.0);
    Eigen::Matrix<double, 4, 2> B;
    B << 1.0, 0.0,
         1.0, 0.0,
         0.0,-1.0,
         0.0,-1.0;
    Eigen::Vector4d z = Eigen::Vector4d::Zero();
    const Eigen::Vector2d x = solver().solve(
        H, f, B, z, Eigen::Vector2d::Zero());
    return x.allFinite() && near(x, Eigen::Vector2d::Zero(), 1e-8)
        && (B*x-z).maxCoeff() <= 1e-8;
}

void active_bound_problem(Eigen::Matrix2d &H,
                          Eigen::Vector2d &f,
                          Eigen::Matrix<double, 1, 2> &B,
                          Eigen::Vector<double, 1> &z,
                          Eigen::Vector2d &x0)
{
    H.setIdentity();
    f << -2.0, 0.5;
    B << 1.0, 0.0;
    z << 0.25;
    x0.setZero();
}

bool truncated_but_feasible_contract()
{
    Eigen::Matrix2d H;
    Eigen::Vector2d f, x0;
    Eigen::Matrix<double, 1, 2> B;
    Eigen::Vector<double, 1> z;
    active_bound_problem(H, f, B, z, x0);
    SolverOptions<double> options;
    options.method = "active set";
    options.stepSizeTolerance = 1e-10;
    options.maxSteps = 1;
    QPSolver<double> qp(options);
    const Eigen::Vector2d x = qp.solve(H, f, B, z, x0);
    const auto result = qp.results();
    return x.allFinite() && (B*x-z).maxCoeff() <= 1e-12
        && !result.converged && result.numberOfSteps == options.maxSteps;
}

bool genuine_convergence_contract()
{
    Eigen::Matrix2d H;
    Eigen::Vector2d f, x0;
    Eigen::Matrix<double, 1, 2> B;
    Eigen::Vector<double, 1> z;
    active_bound_problem(H, f, B, z, x0);
    QPSolver<double> qp = solver();
    const Eigen::Vector2d x = qp.solve(H, f, B, z, x0);
    const Eigen::Vector2d expected(0.25, -0.5);
    return qp.results().converged && near(x, expected, 1e-10);
}

bool convergence_state_does_not_leak()
{
    SolverOptions<double> options;
    options.method = "active set";
    options.stepSizeTolerance = 1e-10;
    options.maxSteps = 2;
    QPSolver<double> qp(options);

    Eigen::Matrix<double, 1, 1> easyH;
    Eigen::Vector<double, 1> easyF, easyX0;
    Eigen::Matrix<double, 0, 1> easyB;
    Eigen::Vector<double, 0> easyZ;
    easyH << 1.0;
    easyF << -1.0;
    easyX0 << 0.0;
    (void)qp.solve(easyH, easyF, easyB, easyZ, easyX0);
    if(!qp.results().converged) return false;

    Eigen::Matrix2d H;
    Eigen::Vector2d f, x0;
    Eigen::Matrix<double, 1, 2> B;
    Eigen::Vector<double, 1> z;
    active_bound_problem(H, f, B, z, x0);
    const Eigen::Vector2d x = qp.solve(H, f, B, z, x0);
    return x.allFinite() && (B*x-z).maxCoeff() <= 1e-12
        && !qp.results().converged
        && qp.results().numberOfSteps == options.maxSteps;
}

bool random_box_qps()
{
    std::mt19937 rng(7);
    std::normal_distribution<double> normal(0.0, 1.0);
    std::uniform_real_distribution<double> uniform(-1.0, 1.0);

    for(int trial = 0; trial < 100; ++trial)
    {
        const int n = 1 + (trial % 5);
        Eigen::MatrixXd M(n, n);
        for(int r = 0; r < n; ++r)
            for(int c = 0; c < n; ++c)
                M(r, c) = normal(rng);
        Eigen::MatrixXd H = M.transpose() * M
                          + 0.1 * Eigen::MatrixXd::Identity(n, n);
        Eigen::VectorXd f(n);
        Eigen::VectorXd lower(n);
        Eigen::VectorXd upper(n);
        Eigen::VectorXd x0(n);
        for(int i = 0; i < n; ++i)
        {
            f(i) = normal(rng);
            lower(i) = -1.5 + 0.25 * uniform(rng);
            upper(i) =  1.5 + 0.25 * uniform(rng);
            x0(i) = lower(i) + (upper(i) - lower(i)) * (0.5 + 0.25 * uniform(rng));
        }

        Eigen::MatrixXd B(2 * n, n);
        Eigen::VectorXd z(2 * n);
        B.topRows(n).setIdentity();
        B.bottomRows(n) = -Eigen::MatrixXd::Identity(n, n);
        z.head(n) = upper;
        z.tail(n) = -lower;

        const Eigen::VectorXd x = solver().solve(H, f, B, z, x0);
        if(!x.allFinite() || !check_kkt(H, f, B, z, x, 1e-7))
        {
            std::cerr << "random QP KKT check failed on trial " << trial
                      << " with x=" << x.transpose() << '\n';
            return false;
        }
    }
    return true;
}

} // namespace

int main()
{
    if(!unconstrained_pd_qp()) return 1;
    if(!active_bound_optimum()) return 2;
    if(!releases_initially_active_constraint()) return 3;
    if(!redundant_boundary_seed_is_well_posed()) return 4;
    if(!truncated_but_feasible_contract()) return 5;
    if(!genuine_convergence_contract()) return 6;
    if(!convergence_state_does_not_leak()) return 7;
    if(!random_box_qps()) return 8;
    return 0;
}
