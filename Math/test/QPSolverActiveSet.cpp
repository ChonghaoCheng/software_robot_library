#include <Math/QPSolver.h>

#include <Eigen/Core>

#include <cmath>

int main()
{
    SolverOptions<double> options;
    options.stepSizeTolerance = 0.05;
    options.maxSteps = 10;
    QPSolver<double> solver(options);

    Eigen::Matrix2d H = Eigen::Matrix2d::Identity();
    Eigen::Vector2d f(-1.0, -1.0);
    Eigen::Matrix<double, 1, 2> B;
    B << 1.0, 0.0;
    Eigen::Vector<double, 1> z;
    z << 0.5;
    Eigen::Vector2d x0(0.5 - 1e-12, 0.0);
    const Eigen::Vector2d blocked = solver.solve(H, f, B, z, x0);
    if(std::abs(blocked.x() - 0.5) > 1e-8 or std::abs(blocked.y() - 1.0) > 1e-8)
        return 1;

    Eigen::Matrix<double, 1, 1> H1;
    H1 << 1.0;
    Eigen::Vector<double, 1> f1;
    f1 << 0.0;
    Eigen::Matrix<double, 1, 1> B1;
    B1 << 1.0;
    Eigen::Vector<double, 1> z1;
    z1 << 0.5;
    Eigen::Vector<double, 1> boundary;
    boundary << 0.5;
    const Eigen::Vector<double, 1> released = solver.solve(H1, f1, B1, z1, boundary);
    return std::abs(released.x()) <= 1e-8 ? 0 : 2;
}
