/**
 * @file    Line.tpp
 * @author  Jon Woolfrey
 * @email   jonathan.woolfrey@gmail.com
 * @date    September 2025
 * @version 1.0
 * @brief   Source file for the Line class.
 *
 * @details A line is defined by a point p0 on the line and a direction vector d. Any point p on the
 *          line satisfies p = p0 + λ d. The point_on_surface() method returns the closest point
 *          on the line to a given reference point.
 *
 * @copyright (c) 2025 Jon Woolfrey
 *
 * @license   OSCL - Free for non-commercial open-source use only.
 *            Commercial use requires a license.
 * 
 * @see https://github.com/Woolfrey/software_robot_library for more information.
 */

#include <cassert>

namespace RobotLibrary { namespace Math {

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                           Constructor                                           //
////////////////////////////////////////////////////////////////////////////////////////////////////
template <unsigned int Dim>
Line<Dim>::Line(const VectorType &direction)
: _direction(direction.normalized())
{
    assert(direction.norm() > 0 && "[ERROR] [LINE] Constructor: Direction vector must be non-zero.");
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                           Get closest point on the line to reference point                     //
////////////////////////////////////////////////////////////////////////////////////////////////////
template <unsigned int Dim>
typename Line<Dim>::VectorType
Line<Dim>::point_on_surface(const VectorType &referencePoint) const
{
    double distance = _direction.dot(referencePoint);                                               // Get distance along the direction vector
    
    return  distance * _direction;                                                                  // Return the point
}

} } // namespace

