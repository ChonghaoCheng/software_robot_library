#include <Math/CondensedMPC.h>

#include <Eigen/Core>
#include <iostream>
#include <vector>

int main()
{
    std::vector<Eigen::MatrixXd> A(3, Eigen::Matrix2d::Identity());
    std::vector<Eigen::MatrixXd> B(3, Eigen::Matrix<double,2,1>::Zero());
    A[0] << 1.0, 0.2, 0.0, 1.0;
    A[1] << 0.9, 0.0, 0.1, 1.1;
    A[2] << 1.0, -0.1, 0.2, 0.8;
    B[0] << 0.1, 0.3;
    B[1] << -0.2, 0.4;
    B[2] << 0.5, 0.1;

    const auto condensed =
        RobotLibrary::Math::condense_time_varying_prediction(A, B);
    const Eigen::Vector2d x0(0.3, -0.4);
    Eigen::Vector3d U;
    U << 0.2, -0.1, 0.5;
    const Eigen::VectorXd stacked =
        condensed.stateTransition * x0 + condensed.inputResponse * U;

    Eigen::Vector2d x = x0;
    for(int k = 0; k < 3; ++k)
    {
        x = A[static_cast<size_t>(k)] * x
          + B[static_cast<size_t>(k)] * U(k);
        const double residual = (stacked.segment<2>(2 * k) - x).norm();
        if(residual > 1e-13)
        {
            std::cerr << "time-varying condensation residual: " << residual << '\n';
            return 1;
        }
    }
    return 0;
}
