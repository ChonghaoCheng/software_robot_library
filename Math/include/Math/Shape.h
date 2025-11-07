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
         * @brief Get a point on the surface / circumference of the shape.
         * @note Not always the closest point.
         * @param referencePoint An external reference point used for computation.
         * @return What you asked for.
         */
        virtual
        VectorType
        point_on_surface(const VectorType &referencePoint) const = 0;

        virtual
        Eigen::Matrix2d
        get_shape_matrix() const = 0;

};

// Convenience aliases
using Shape2D = Shape<2>;
using Shape3D = Shape<3>;

} } // namespace RobotLibrary::Math

#endif

