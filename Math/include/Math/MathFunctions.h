/**
 * @file    MathFunctions.h
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

#ifndef MATH_FUNCTIONS_H
#define MATH_FUNCTIONS_H

#include <Math/DataStructures.h>

#include <Eigen/Dense>                                                                               // Eigen::Vector, Eigen::Matrix etc
#include <Eigen/Geometry>                                                                            // Eigen::Quaterniond
#include <iostream>
#include <vector>

namespace RobotLibrary { namespace Math {

/**
 * It's obvious what this function does.
 * @param A a square matrix
 * @return True if positive-definite, false otherwise
 */
bool is_positive_definite(const Eigen::MatrixXd &A,
                          std::string &info);
       
/**
 * Decompose a matrix A = Q*R where Q is an orthogonal matrix, and R is upper-triangular.
 * @param A The matrix to be decomposed.
 * @param tolerance The rounding error on a singularity
 * @return A QRDecomposition data structure
 */
RobotLibrary::Math::QRDecomposition
schwarz_rutishauser(const Eigen::MatrixXd &A,
                    const double tolerance = 1e-04);

/**
 * Solve a system of equations y = L*x, where L is a lower-triangular matrix.
 * @param L A lower-triangular matrix.
 * @param U A vector of known values.
 * @param tolerance For singularities.
 * @return A solution for x.
 */
Eigen::MatrixXd
forward_substitution(const Eigen::MatrixXd &L,
                     const Eigen::MatrixXd &Y,
                     const double tolerance = 1e-04);

/**
 * Solve a system of equations Y = U*X, where U is an upper-triangular matrix.
 * @param U An upper-triangular matrix.
 * @param Y A tensor (vector or matrix) of known values.
 * @param tolerance For handling singularities.
 * @return A solution for X.
 */
Eigen::MatrixXd
backward_substitution(const Eigen::MatrixXd &U,
                      const Eigen::MatrixXd &Y,         
                      const double tolerance = 1e-04);

/**
 * This function solves for the derivatives at each point of a cubic spline such that there is continuity.
 * @param y The dependent variable for the spline.
 * @param x The independent variable for the spline.
 * @param firstDerivative The value for dy/dx at the very first point.
 * @param finalDerivative The value for dy/dx at the final point.
 * @return An array containing the derivatives dy/dx for each point x.
 */
std::vector<double>
solve_cubic_spline_derivatives(const std::vector<double> &y,
                               const std::vector<double> &x,
                               const double &firstDerivative = 0,
                               const double &finalDerivative = 0);

/**
 * @brief Ensure that an angle is between -3.141592.. and 3.141592...
 * @param angle The value to be capped.
 * @return The angle re-mapped.
 */
double
wrap_to_pi(const double &angle);

/**
 * @brief Logarithmic map of SO(3): a rotation matrix to a rotation vector (angle * axis).
 * @param R A 3x3 rotation matrix.
 * @return The rotation vector in R^3.
 */
Eigen::Vector3d
so3_logarithm(const Eigen::Matrix3d &R);

/**
 * @brief Logarithmic map of SE(3): a homogeneous transform to a body twist.
 * @param T A 4x4 homogeneous transform.
 * @return The body twist [v; w] in R^6 with the [vx vy vz wx wy wz] ordering.
 */
Eigen::Matrix<double,6,1>
se3_logarithm(const Eigen::Matrix4d &T);

/**
 * @brief Lie-algebra adjoint ad_xi for an SE(3) twist ordered [v; w].
 * @details ad_xi * eta is the Lie bracket [xi, eta].
 */
Eigen::Matrix<double,6,6>
se3_adjoint_matrix(const Eigen::Matrix<double,6,1> &twist);

/**
 * @brief Right Jacobian of SE(3), mapping log-coordinate rate to body twist.
 * @details xi_body = J_r(e) * e_dot for T = Exp(e).
 */
Eigen::Matrix<double,6,6>
se3_right_jacobian(const Eigen::Matrix<double,6,1> &twist);

/**
 * @brief Inverse right Jacobian of SE(3): e_dot = J_r(e)^-1 * xi_body.
 * @throws std::runtime_error if the Jacobian is numerically singular.
 */
Eigen::Matrix<double,6,6>
se3_right_jacobian_inverse(const Eigen::Matrix<double,6,1> &twist);

/**
 * @brief Inverse of a homogeneous SE(3) transform (uses the block structure, not a generic inverse).
 * @param T A 4x4 homogeneous transform.
 * @return The inverse transform.
 */
Eigen::Matrix4d
se3_inverse(const Eigen::Matrix4d &T);

/**
 * @brief Exponential map of SO(3): a rotation vector (angle * axis) to a rotation matrix.
 *        Inverse of so3_logarithm(). Uses Rodrigues' formula with small-angle Taylor expansion.
 * @param vector The rotation vector in R^3.
 * @return A 3x3 rotation matrix.
 */
Eigen::Matrix3d
so3_exponential(const Eigen::Vector3d &vector);

/**
 * @brief Exponential map of SE(3): a body twist to a homogeneous transform.
 *        Inverse of se3_logarithm(). The twist uses the [v; w] (linear first) ordering.
 * @param twist The body twist [v; w] in R^6.
 * @return A 4x4 homogeneous transform.
 */
Eigen::Matrix4d
se3_exponential(const Eigen::Matrix<double,6,1> &twist);

/**
 * @brief Left Jacobian of SO(3). Maps a rotation-vector rate to the SPATIAL angular velocity:
 *        w_spatial = J_l(r) * r_dot. Equivalently the SE(3) "V" matrix used by se3_exponential().
 * @param vector The rotation vector r in R^3.
 * @return A 3x3 matrix.
 */
Eigen::Matrix3d
so3_left_jacobian(const Eigen::Vector3d &vector);

/**
 * @brief Right Jacobian of SO(3). Maps a rotation-vector rate to the BODY angular velocity:
 *        w_body = J_r(r) * r_dot. Note J_r(r) = J_l(-r).
 * @param vector The rotation vector r in R^3.
 * @return A 3x3 matrix.
 */
Eigen::Matrix3d
so3_right_jacobian(const Eigen::Vector3d &vector);

/**
 * @brief Inverse of the left Jacobian of SO(3): r_dot = J_l^{-1}(r) * w_spatial.
 *        This is exactly the "Vinv" matrix used inside se3_logarithm().
 * @param vector The rotation vector r in R^3.
 * @return A 3x3 matrix.
 */
Eigen::Matrix3d
so3_left_jacobian_inverse(const Eigen::Vector3d &vector);

/**
 * @brief Inverse of the right Jacobian of SO(3): r_dot = J_r^{-1}(r) * w_body.
 * @param vector The rotation vector r in R^3.
 * @return A 3x3 matrix.
 */
Eigen::Matrix3d
so3_right_jacobian_inverse(const Eigen::Vector3d &vector);

/**
 * @brief Adjoint of an SE(3) transform for the [v; w] (linear first) twist ordering.
 *        Maps a body twist to a spatial twist: xi_spatial = Adjoint(T) * xi_body, i.e.
 *        v_spatial = R*v_body + p x (R*w_body),  w_spatial = R*w_body.
 * @param T A 4x4 homogeneous transform.
 * @return A 6x6 adjoint matrix.
 */
Eigen::Matrix<double,6,6>
adjoint(const Eigen::Matrix4d &T);

/**
 * @brief Convert a quaternion to the shortest equivalent rotation vector (angle * axis).
 *        Robust to non-normalised / non-finite inputs, and enforces the w >= 0 hemisphere.
 * @param quaternion The orientation to convert.
 * @return The rotation vector in R^3.
 */
Eigen::Vector3d
quaternion_to_rotation_vector(const Eigen::Quaterniond &quaternion);

/**
 * @brief Orientation error from a current to a desired orientation, as a rotation vector
 *        expressed in the reference (world) frame: e = log( q_desired * q_current^{-1} ).
 * @param current The current orientation.
 * @param desired The desired orientation.
 * @return The orientation error rotation vector in R^3.
 */
Eigen::Vector3d
quaternion_orientation_error(const Eigen::Quaterniond &current,
                             const Eigen::Quaterniond &desired);

} }

#endif
