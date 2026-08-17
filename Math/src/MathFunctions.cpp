/**
 * @file    MathFunctions.cpp
 * @author  Jon Woolfrey
 * @email   jonathan.woolfrey@gmail.com
 * @date    July 2025
 * @version 1.1
 * @brief   Useful math functions for robot kinematics & control.
 * 
 * @details This header files contains forward declarations for useful math functions that are not
 *          offered by the Eigen library.
 * 
 * @copyright (c) 2025 Jon Woolfrey
 *
 * @license   OSCL - Free for non-commercial open-source use only.
 *            Commercial use requires a license.
 * 
 * @see https://github.com/Woolfrey/software_robot_library for more information.
 */
 
#include <Math/MathFunctions.h>
#include <Math/SkewSymmetric.h>

#include <cmath>

namespace RobotLibrary { namespace Math {

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                      POSITIVE DEFINITE?                                        //
////////////////////////////////////////////////////////////////////////////////////////////////////
bool is_positive_definite(const Eigen::MatrixXd &A,
                          std::string &message)
{
    if (A.rows() != A.cols())
    {
        message = "Matrix is not square.";
        return false;
    }

    double symmetryErrorNorm = (A - A.transpose()).norm();
    if (symmetryErrorNorm > 1e-4)
    {
        message = "Matrix is not symmetric; ||A - A'|| = " + std::to_string(symmetryErrorNorm) + " > 1e-4.\n";
        return false;
    }

    Eigen::LLT<Eigen::MatrixXd> llt(A);
    if (llt.info() == Eigen::Success)
    {
        return true;
    }
    else
    {
        message = "Cholesky decomposition failed.";
        return false;
    }
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                     QR DECOMPOSITION                                           //
////////////////////////////////////////////////////////////////////////////////////////////////////
RobotLibrary::Math::QRDecomposition
schwarz_rutishauser(const Eigen::MatrixXd &A,
                    const double tolerance)
{
     unsigned int m = A.rows();
     unsigned int n = A.cols();
     
     if(m < n)
     {
          throw std::invalid_argument("[ERROR] qr_decomposition() "
                                      "Matrix A has " + std::to_string(A.rows()) +  " rows which is "
                                      "less than its " + std::to_string( A.cols()) +  " columns. "
                                      "Cannot solve the QR decomposition.");
     }
     else
     {
          // Schwarz-Rutishauser Algorithm.
          // A full decomposition of an mxn matrix A (m > n) is:
          //    [ Qr Qn ][ R ]  = A
          //             [ 0 ]
          //
          // where: - Qr is mxn,
          //        - Qn is mx(m-n)
          //        - R  is nxn
          //
          // The null space of A is obtained with N = Qn*Qn'.
          // This algorithm returns only Qr and R for efficiency.
          
          RobotLibrary::Math::QRDecomposition decomp;
          decomp.Q = A;
          decomp.R.resize(n,n); decomp.R.setZero();
          
          for(int j = 0; j < n; j++)
          {
               for(int i = 0; i < j; i++)
               {
                    decomp.R(i,j)   = decomp.Q.col(i).dot(decomp.Q.col(j));                         // Project the columns
                    decomp.Q.col(j) = decomp.Q.col(j) - decomp.R(i,j)*decomp.Q.col(i);
               }
               
               decomp.R(j,j) = decomp.Q.col(j).norm();
               
               if(abs(decomp.R(j,j)) > tolerance) decomp.Q.col(j) /= decomp.R(j,j);
               else                               decomp.Q.col(j).setZero();                        // Singular
          }
          
          return decomp;
     }
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                    FORWARD SUBSTITUTION                                        //
////////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::MatrixXd
forward_substitution(const Eigen::MatrixXd &Y,
                     const Eigen::MatrixXd &L,
                     const double tolerance)
{
     unsigned int m = Y.rows();
     unsigned int n = L.cols();
     unsigned int o = Y.cols();
     
     if(L.rows() != m)
     {
          throw std::invalid_argument("[ERROR] forward_subsitution(): "
                                      "Dimensions of arguments do not match. "
                                      "The vector input had " + std::to_string(m) + " elements, "
                                      "but the matrix input had " + std::to_string(L.rows()) + " rows.");
     }
     else if(L.rows() != L.cols())
     {
          throw std::invalid_argument("[ERROR] forward_substiution(): "
                                      "Expected a square matrix but it was " + std::to_string(L.rows()) +
                                      "x" + std::to_string(L.cols()) + ".");
     }
     
     Eigen::MatrixXd X(m,o);                                                                        // Value to be returned
     
     for(int i = 0; i < o; i++)
     {
          for(int j = 0; j < m; j++)
          {
               double sum = 0.0;
               
               for(int k = 0; k < j; k++) sum += L(j,k)*X(k,i);
               
               if(L(j,j) >= tolerance) X(j,i) = (Y(j,i) - sum)/L(j,j);
               else                    X(j,i) = 0;
          }
     }
     
     return X;
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                    BACKWARD SUBSTITUTION                                       //
////////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::MatrixXd
backward_substitution(const Eigen::MatrixXd &U,
                      const Eigen::MatrixXd &Y,
                      const double tolerance)
{
     unsigned int m = Y.rows();
     unsigned int n = U.cols();
     unsigned int o = Y.cols();
       
     if(U.rows() != m)
     {
          throw std::invalid_argument("[ERROR] forward_subsitution(): "
                                      "Dimensions of arguments do not match. "
                                      "The vector input had " + std::to_string(m) + " rows, "
                                      "but the matrix input had " + std::to_string(U.rows()) + " rows.");
     }
     else if(U.rows() != U.cols())
     {
          throw std::invalid_argument("[ERROR] forward_substiution(): "
                                      "Expected a square matrix but it was " + std::to_string(U.rows()) +
                                      "x" + std::to_string(U.cols()) + ".");
     }
     
     Eigen::MatrixXd X(n,o);                                                                        // Value to be turned
     
     for(int i = 0; i < o; i++)                                                                     // For every column of Y
     {
          for(int j = m-1; j >= 0; j--)
          {
               double sum = 0.0;
               
               for(int k = j + 1; k < m; k++) sum += U(j,k)*X(k,i);
               
               if(U(j,j) >= tolerance) X(j,i) = (Y(j,i)-sum)/U(j,j);
               else                    X(j,i) = 0.0;
          }
     }
     
     return X;
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                        Fit the derivatives for cubic spline interpolation                      //
////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double>
solve_cubic_spline_derivatives(const std::vector<double> &y,
                               const std::vector<double> &x,
                               const double &firstDerivative,
                               const double &finalDerivative)
{
    unsigned int n = y.size();
    
    if(n < 3)
    {
        throw std::invalid_argument("[ERROR] solve_cubic_spline derivatives(): "
                                    "A minimum number of 3 points is required to define a spline.");
    }
    else if(y.size() != x.size())
    {
        throw std::invalid_argument("[ERROR] solve_cubic_spline_derivatives(): "
                                    "Dimensions of arguments do not match. The y vector had " +
                                    std::to_string(y.size()) + " elements, and the x vector had " +
                                    std::to_string(x.size()) + " elements.");
    }
    
    using namespace Eigen;
    
    // Derivatives at waypoints are related to positions via the relationship:
    // A*dy = B*y --> dy = A^-1*B*y
    
    Eigen::MatrixXd A(n,n); A.setZero();
    
    Eigen::MatrixXd B = A;
    
    A(0,0) = 1;                                                                                     // First value
    
    // Assign intermediate values
    for(int i = 1; i < n-1; i++)
    {
        double dx1 = x[i]   - x[i-1];
        double dx2 = x[i+1] - x[i];

        if(dx1 == 0)
        {
            throw std::logic_error("[ERROR] fit_cubic_spline_derivatives(): "
                                   "Independent variable " + std::to_string(i) + " is the same as "
                                   "independent variable " + std::to_string(i+1) + " ("
                                   + std::to_string(x[i-1]) + " == " + std::to_string(x[i]) + ").");
        }

        A(i,i-1) = 1/dx1;
        A(i,i)   = 2*(1/dx1 + 1/dx2);
        A(i,i+1) = 1/dx2;

        B(i,i-1) = -3/(dx1*dx1);
        B(i,i)   =  3*(1/(dx1*dx1) - 1/(dx2*dx2));
        B(i,i+1) =  3/(dx2*dx2);
    }

    A(n-1,n-1) = 1;
    
    Eigen::VectorXd points(y.size());
    for(int i = 0; i < y.size(); i++) points(i) = y[i];
    
    Eigen::VectorXd derivatives = A.partialPivLu().solve(B*points);
    
    std::vector<double> temp(derivatives.size());
    
    for(int i = 0; i < derivatives.size(); i++) temp[i] = derivatives[i];
    
    return temp;                                  
}

  ///////////////////////////////////////////////////////////////////////////////////////////////////
 //                          Map an angle to the interveral (-pi, pi]                             //
///////////////////////////////////////////////////////////////////////////////////////////////////
double
wrap_to_pi(const double &angle)
{
    double newAngle = fmod(angle + M_PI, 2 * M_PI);

    if(newAngle < 0) newAngle += 2 * M_PI;

    return newAngle - M_PI;
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                            Logarithmic map of SO(3): R -> rotation vector                       //
////////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Vector3d
so3_logarithm(const Eigen::Matrix3d &R)
{
    Eigen::Quaterniond q(R);
    q.normalize();
    if(q.w() < 0.0) q.coeffs() *= -1.0;                                                             // Shortest path

    const Eigen::Vector3d vector = q.vec();
    const double vectorNorm = vector.norm();
    if(not std::isfinite(vectorNorm) or vectorNorm < 1e-10) return Eigen::Vector3d::Zero();

    const double angle = 2.0 * std::atan2(vectorNorm, q.w());
    return angle * vector / vectorNorm;
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                            Logarithmic map of SE(3): T -> body twist                            //
////////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Matrix<double,6,1>
se3_logarithm(const Eigen::Matrix4d &T)
{
    Eigen::Matrix<double,6,1> twist = Eigen::Matrix<double,6,1>::Zero();
    const Eigen::Matrix3d R = T.block<3,3>(0,0);
    const Eigen::Vector3d p = T.block<3,1>(0,3);
    const Eigen::Vector3d w = so3_logarithm(R);

    // The translation maps through the inverse left Jacobian: v = J_l^{-1}(w) * p.
    twist.head<3>() = so3_left_jacobian_inverse(w) * p;
    twist.tail<3>() = w;
    return twist;
}

Eigen::Matrix<double,6,6>
se3_adjoint_matrix(const Eigen::Matrix<double,6,1> &twist)
{
    const Eigen::Matrix3d V = SkewSymmetric(twist.head<3>()).as_matrix();
    const Eigen::Matrix3d W = SkewSymmetric(twist.tail<3>()).as_matrix();

    Eigen::Matrix<double,6,6> ad = Eigen::Matrix<double,6,6>::Zero();
    ad.block<3,3>(0,0) = W;
    ad.block<3,3>(0,3) = V;
    ad.block<3,3>(3,3) = W;
    return ad;
}

Eigen::Matrix<double,6,6>
se3_right_jacobian(const Eigen::Matrix<double,6,1> &twist)
{
    using Matrix6d = Eigen::Matrix<double,6,6>;
    const Matrix6d negativeAd = -se3_adjoint_matrix(twist);

    // J_r(e) = sum_{n=0}^infinity (-ad_e)^n / (n+1)!.
    // The logarithm uses the shortest SO(3) branch, so a short fixed series is
    // accurate throughout the controller's nonsingular operating region.
    Matrix6d jacobian = Matrix6d::Identity();
    Matrix6d term = Matrix6d::Identity();
    for(int order = 1; order <= 24; ++order)
    {
        term = term * negativeAd / static_cast<double>(order + 1);
        jacobian += term;
        if(term.norm() < 1e-15) break;
    }
    return jacobian;
}

Eigen::Matrix<double,6,6>
se3_right_jacobian_inverse(const Eigen::Matrix<double,6,1> &twist)
{
    using Matrix6d = Eigen::Matrix<double,6,6>;
    const Matrix6d jacobian = se3_right_jacobian(twist);
    const Eigen::FullPivLU<Matrix6d> decomposition(jacobian);
    if(not decomposition.isInvertible())
    {
        throw std::runtime_error(
            "[ERROR] [MATH FUNCTIONS] se3_right_jacobian_inverse(): Jacobian is singular.");
    }
    return decomposition.solve(Matrix6d::Identity());
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                              Inverse of a homogeneous SE(3) transform                           //
////////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Matrix4d
se3_inverse(const Eigen::Matrix4d &T)
{
    Eigen::Matrix4d inverse = Eigen::Matrix4d::Identity();
    const Eigen::Matrix3d R = T.block<3,3>(0,0);
    const Eigen::Vector3d p = T.block<3,1>(0,3);
    inverse.block<3,3>(0,0) = R.transpose();
    inverse.block<3,1>(0,3) = -R.transpose() * p;
    return inverse;
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                            Exponential map of SO(3): rotation vector -> R                       //
////////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Matrix3d
so3_exponential(const Eigen::Vector3d &vector)
{
    const double theta  = vector.norm();
    const Eigen::Matrix3d W = SkewSymmetric(vector).as_matrix();

    // Rodrigues: R = I + a*W + b*W^2, with a = sin(t)/t, b = (1 - cos(t))/t^2.
    // Use Taylor series near 0 to avoid the 0/0 singularity (and keep R orthonormal to O(t^4)).
    double a, b;
    if(theta < 1e-6)
    {
        const double t2 = theta * theta;
        a = 1.0 - t2 / 6.0;                                                                          // sin(t)/t
        b = 0.5 - t2 / 24.0;                                                                         // (1 - cos t)/t^2
    }
    else
    {
        a = std::sin(theta) / theta;
        b = (1.0 - std::cos(theta)) / (theta * theta);
    }

    return Eigen::Matrix3d::Identity() + a * W + b * W * W;
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                              Left Jacobian of SO(3)  (and inverse)                              //
////////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Matrix3d
so3_left_jacobian(const Eigen::Vector3d &vector)
{
    const double theta = vector.norm();
    const Eigen::Matrix3d W = SkewSymmetric(vector).as_matrix();

    // J_l = I + c1*W + c2*W^2,  c1 = (1 - cos t)/t^2,  c2 = (t - sin t)/t^3
    double c1, c2;
    if(theta < 1e-6)
    {
        const double t2 = theta * theta;
        c1 = 0.5 - t2 / 24.0;
        c2 = 1.0 / 6.0 - t2 / 120.0;
    }
    else
    {
        const double t2 = theta * theta;
        c1 = (1.0 - std::cos(theta)) / t2;
        c2 = (theta - std::sin(theta)) / (t2 * theta);
    }

    return Eigen::Matrix3d::Identity() + c1 * W + c2 * W * W;
}

Eigen::Matrix3d
so3_right_jacobian(const Eigen::Vector3d &vector)
{
    return so3_left_jacobian(-vector);                                                               // J_r(r) = J_l(-r)
}

Eigen::Matrix3d
so3_left_jacobian_inverse(const Eigen::Vector3d &vector)
{
    const double theta = vector.norm();
    const Eigen::Matrix3d W = SkewSymmetric(vector).as_matrix();

    // J_l^{-1} = I - 0.5*W + c*W^2,  c = 1/t^2 - (1 + cos t)/(2 t sin t)  (-> 1/12 as t -> 0)
    // The (1 + cos t)/sin t form cancels catastrophically near t = pi (both -> 0), so use the
    // equivalent and well-conditioned half-angle form (1 + cos t)/sin t = cot(t/2) instead.
    double c;
    if(theta < 1e-6)
    {
        c = 1.0 / 12.0 + theta * theta / 720.0;
    }
    else
    {
        const double half = 0.5 * theta;
        c = 1.0 / (theta * theta) - std::cos(half) / (2.0 * theta * std::sin(half));
    }

    return Eigen::Matrix3d::Identity() - 0.5 * W + c * W * W;
}

Eigen::Matrix3d
so3_right_jacobian_inverse(const Eigen::Vector3d &vector)
{
    return so3_left_jacobian_inverse(-vector);                                                       // J_r^{-1}(r) = J_l^{-1}(-r)
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                            Exponential map of SE(3): body twist -> T                            //
////////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Matrix4d
se3_exponential(const Eigen::Matrix<double,6,1> &twist)
{
    const Eigen::Vector3d v = twist.head<3>();
    const Eigen::Vector3d w = twist.tail<3>();

    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T.block<3,3>(0,0) = so3_exponential(w);
    T.block<3,1>(0,3) = so3_left_jacobian(w) * v;                                                    // p = V * v, with V = J_l (matches se3_logarithm's Vinv = J_l^{-1})
    return T;
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                  Adjoint of an SE(3) transform                                  //
////////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Matrix<double,6,6>
adjoint(const Eigen::Matrix4d &T)
{
    const Eigen::Matrix3d R = T.block<3,3>(0,0);
    const Eigen::Vector3d p = T.block<3,1>(0,3);

    Eigen::Matrix<double,6,6> Ad = Eigen::Matrix<double,6,6>::Zero();
    Ad.block<3,3>(0,0) = R;
    Ad.block<3,3>(0,3) = SkewSymmetric(p).as_matrix() * R;                                           // [v;w] ordering: top-right block is [p]x * R
    Ad.block<3,3>(3,3) = R;
    return Ad;
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                              Quaternion -> shortest rotation vector                             //
////////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Vector3d
quaternion_to_rotation_vector(const Eigen::Quaterniond &quaternion)
{
    Eigen::Quaterniond q = quaternion;
    q.normalize();

    if(not std::isfinite(q.w()) or not q.vec().allFinite()) return Eigen::Vector3d::Zero();

    if(q.w() < 0.0) q.coeffs() *= -1.0;                                                              // Shortest path (w >= 0 hemisphere)

    const Eigen::Vector3d v = q.vec();
    const double vNorm = v.norm();
    if(vNorm < 1e-8) return Eigen::Vector3d::Zero();

    const double angle = 2.0 * std::atan2(vNorm, q.w());
    return angle * v / vNorm;
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                       Orientation error (world frame) as a rotation vector                      //
////////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Vector3d
quaternion_orientation_error(const Eigen::Quaterniond &current,
                             const Eigen::Quaterniond &desired)
{
    return quaternion_to_rotation_vector(desired * current.conjugate());                             // log( q_d * q_c^{-1} )
}

} }
