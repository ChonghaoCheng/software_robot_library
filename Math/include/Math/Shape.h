/**
 * @file    Shape.h
 * @author  Jon Woolfrey
 * @email   jonathan.woolfrey@gmail.com
 * @date    September 2025
 * @version 1.0
 * @brief   A class used for representing shapes (duh!)
 *
 * @details This abstract class provides a common interface for all shapes.
 *
 * @copyright (c) 2025 Jon Woolfrey
 *
 * @license   OSCL - Free for non-commercial open-source use only.
 *            Commercial use requires a license.
 * 
 * @see https://github.com/Woolfrey/software_robot_library for more information.
 */

#ifndef SHAPE_H
#define SHAPE_H

#include <Eigen/Dense>

namespace RobotLibrary { namespace Math {

template <unsigned int Dim>
class Shape
{
    public:
    
        using VectorType = Eigen::Vector<double, Dim>;                                              // This is for brevity

        /**
         * @brief Destructor.
         */
        virtual
        ~Shape() = default;

        /**
         * @brief Computes the closest distance to the surface of the shape.
         * @param point A point in n-dimensional space.
         * @return A scalar value.
         * @note This is a virtual method and must be defined in any child class.
         */
        virtual
        double
        distance_to_surface(const VectorType& point) const = 0;

        /**
         * @brief Computes the squared distance to the closest point on the surface of the shape.
         * @param point A point in n-dimensional space.
         * @return A positive scalar.
         * @note This virtual method should be defined in any child class.
         */
        virtual
        double
        distance_to_surface_squared(const VectorType& point) const = 0;
};

// Convenience aliases
using Shape2D = Shape<2>;
using Shape3D = Shape<3>;

} } // namespace RobotLibrary::Math

#endif

