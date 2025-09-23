/**
 * @file    Line.h
 * @author  Jon Woolfrey
 * @email   jonathan.woolfrey@gmail.com
 * @date    September 2025
 * @version 1.0
 * @brief   A class that represents a line in n-dimensional space.
 *
 * @details A line is defined by a point p0 on the line, and a direction vector d. Any point p on the
 *          line satisfies p = p0 + λ d for some scalar λ.
 *
 * @copyright (c) 2025 Jon Woolfrey
 *
 * @license   OSCL - Free for non-commercial open-source use only.
 *            Commercial use requires a license.
 * 
 * @see https://github.com/Woolfrey/software_robot_library for more information.
 */

#ifndef LINE_H
#define LINE_H

#include <Math/Shape.h>                                                                             // Base class

#include <Eigen/Dense>                                                                              // For vector algebra

namespace RobotLibrary { namespace Math {

/**
 * @brief A class for representing n-dimensional lines.
 * @note  A line is defined by a point and a direction vector.
 */
template <unsigned int Dim>
class Line : public Shape<Dim>
{
    public:

        using VectorType = Eigen::Matrix<double, Dim, 1>;

        /**
         * @brief Constructor.
         * @param point A point lying on the line.
         * @param direction A direction vector (will be normalized).
         */
        Line(const VectorType &direction);

        /**
         * @brief Get the closest point on the line to a given reference point.
         * @param referencePoint External point.
         * @return Closest point on the line to referencePoint.
         */
        VectorType
        point_on_surface(const VectorType &referencePoint) const override;

    private:

        VectorType _direction;                                                                      ///< Normalized direction vector.
};

// Convenience aliases
using Line2D = Line<2>;
using Line3D = Line<3>;

} } // namespace

#include <Math/Line.tpp>

#endif
