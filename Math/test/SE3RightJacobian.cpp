#include <Math/MathFunctions.h>

#include <Eigen/Core>
#include <iostream>

int main()
{
    using Vector6d = Eigen::Matrix<double,6,1>;
    using Matrix6d = Eigen::Matrix<double,6,6>;
    using RobotLibrary::Math::adjoint;
    using RobotLibrary::Math::se3_exponential;
    using RobotLibrary::Math::se3_inverse;
    using RobotLibrary::Math::se3_logarithm;
    using RobotLibrary::Math::se3_right_jacobian;
    using RobotLibrary::Math::se3_right_jacobian_inverse;

    const Vector6d zero = Vector6d::Zero();
    if((se3_right_jacobian(zero) - Matrix6d::Identity()).norm() > 1e-14) return 1;
    if((se3_right_jacobian_inverse(zero) - Matrix6d::Identity()).norm() > 1e-14) return 2;

    Vector6d error;
    error << 0.18, -0.12, 0.07, 0.35, -0.22, 0.28;
    Vector6d currentBodyTwist;
    currentBodyTwist << -0.08, 0.04, 0.03, 0.06, -0.02, 0.05;
    Vector6d referenceTangent;
    referenceTangent << 0.24, -0.10, 0.05, 0.12, 0.18, -0.08;
    const double progressRate = 0.31;
    const double dt = 1e-7;

    const Eigen::Matrix4d relative = se3_exponential(error);
    const Eigen::Matrix4d nextRelative =
        se3_exponential(-referenceTangent * progressRate * dt)
        * relative
        * se3_exponential(currentBodyTwist * dt);
    const Vector6d finiteDifference = (se3_logarithm(nextRelative) - error) / dt;
    const Vector6d predicted = se3_right_jacobian_inverse(error)
        * (currentBodyTwist
           - adjoint(se3_inverse(relative)) * referenceTangent * progressRate);

    const double residual = (finiteDifference - predicted).norm();
    if(residual > 2e-7)
    {
        std::cerr << "SE(3) right-Jacobian dynamics residual: " << residual << '\n';
        return 3;
    }
    return 0;
}
