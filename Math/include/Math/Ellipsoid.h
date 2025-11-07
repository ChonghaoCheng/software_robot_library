/**
 * @file    Ellipsoid.h
 * @author  Jon Woolfrey
 * @email   jonathan.woolfrey@gmail.com
 * @date    September 2025
 * @version 1.0
 * @brief   A class that represents an ellipsoid in n-dimensional space.
 *
 * @details An ellipsoid is defined by its centerpoint c, and matrix A such that, for any point p on
 *          it surface it satisfies (p - c)^T A^{-1} (p - c) = 1. The center is always assumed to be
 *          at zero.
 *
 * @copyright (c) 2025 Jon Woolfrey
 *
 * @license   OSCL - Free for non-commercial open-source use only.
 *            Commercial use requires a license.
 * 
 * @see https://github.com/Woolfrey/software_robot_library for more information.
 */

#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include <Math/Shape.h>                                                                             // Base class

#include <Eigen/Dense>                                                                              // Allows Cholesky decomposition

namespace RobotLibrary { namespace Math {

/**
 * @brief A class for representing n-dimensional ellipsoids. Given the center c, and shape matrix A,
 *        an ellipsoid satisfies (p - c)^T * A^-1 * (p - c) = 1 for any point p on its surface.
 * @Note  The center is assumed to be c = 0.
 */
template <unsigned int Dim>
class Ellipsoid : public Shape<Dim>
{
    public:
    
        /**
         * @brief Constructor using a positive-definite matrix.
         * @param center The position of the center.
         * @param shapeMatrix A positive-definite matrix.
         */
        Ellipsoid(const Eigen::Matrix<double,Dim,Dim> &shapeMatrix);
                         
        /**
         * @brief Get a point on the circumference on the ray to the center.
         * @note This method overrides the one in the base class.
         * @param referencePoint An external reference point used for computation.
         * @return What you asked for.
         */
        Eigen::Vector<double,Dim>
        point_on_surface(const Eigen::Vector<double,Dim> &referencePoint) const override;

        Eigen::Matrix<double, Dim, Dim>
        get_shape_matrix() const override;

    private:
        
        Eigen::LLT<Eigen::Matrix<double, Dim, Dim>> _LLT;                                           ///< Cholesky decomposition of the shape matrix.
        
        Eigen::Matrix<double, Dim, Dim> _shapeMatrix;                                               ///< A positive definite matrix describing its shape
};

// Convenience aliases
using Ellipsoid2D = Ellipsoid<2>;
using Ellipsoid3D = Ellipsoid<3>;

} } // namespace

#include <Math/Ellipsoid.tpp>

#endif
