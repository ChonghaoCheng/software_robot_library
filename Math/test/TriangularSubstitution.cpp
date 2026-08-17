#include <Math/MathFunctions.h>

#include <Eigen/Core>

#include <iostream>

int main()
{
    Eigen::Matrix3d upper;
    upper << 2.0, -1.0, 3.0,
             0.0,  4.0, 2.0,
             0.0,  0.0, 5.0;
    Eigen::Matrix<double,3,2> expected;
    expected << 1.0, -2.0,
                3.0,  0.5,
               -1.0,  4.0;
    const Eigen::Matrix<double,3,2> rightHandSide = upper * expected;
    const Eigen::MatrixXd actual =
        RobotLibrary::Math::backward_substitution(upper, rightHandSide, 1e-12);

    if((actual - expected).norm() > 1e-12
       or (upper * actual - rightHandSide).norm() > 1e-12)
    {
        std::cerr << "Backward substitution failed for multiple right-hand sides.\n";
        return 1;
    }
    return 0;
}
