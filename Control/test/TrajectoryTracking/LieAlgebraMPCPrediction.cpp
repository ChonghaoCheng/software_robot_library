#include "LieAlgebraMPCPrediction.h"
#include <Math/MathFunctions.h>

#include <Eigen/Core>
#include <iostream>

int main()
{
    using RobotLibrary::Control::LieAlgebraMPCVector;
    using RobotLibrary::Control::endpoint_twist_in_base;
    using RobotLibrary::Control::endpoint_twist_in_body;
    using RobotLibrary::Control::lie_algebra_mpc_stage;
    using RobotLibrary::Math::se3_exponential;
    using RobotLibrary::Math::se3_logarithm;

    LieAlgebraMPCVector reference;
    reference << 0.12, -0.04, 0.08, 0.16, -0.11, 0.07;
    LieAlgebraMPCVector errorDirection;
    errorDirection << -0.3, 0.2, 0.1, 0.12, -0.08, 0.15;
    LieAlgebraMPCVector correctionDirection;
    correctionDirection << 0.1, -0.15, 0.05, -0.06, 0.09, 0.04;

    const double dt = 1e-5;
    const double epsilon = 1e-5;
    const auto [A, B] = lie_algebra_mpc_stage(reference, dt);
    const LieAlgebraMPCVector predicted =
        A * (epsilon * errorDirection) + B * (epsilon * correctionDirection);

    const Eigen::Matrix4d exact =
        se3_exponential(-dt * reference)
        * se3_exponential(epsilon * errorDirection)
        * se3_exponential(dt * (reference + epsilon * correctionDirection));
    const LieAlgebraMPCVector actual = se3_logarithm(exact);
    const double residual = (actual - predicted).norm();
    if(residual > 1e-9)
    {
        std::cerr << "Lie-algebra stage residual: " << residual << '\n';
        return 1;
    }

    Eigen::AngleAxisd rotation(0.73, Eigen::Vector3d(0.2, -0.4, 0.8).normalized());
    LieAlgebraMPCVector baseTwist;
    baseTwist << 0.2, -0.1, 0.3, 0.4, 0.1, -0.2;
    const LieAlgebraMPCVector roundTrip = endpoint_twist_in_base(
        rotation.toRotationMatrix(),
        endpoint_twist_in_body(rotation.toRotationMatrix(), baseTwist));
    if((roundTrip - baseTwist).norm() > 1e-14) return 2;

    return 0;
}
