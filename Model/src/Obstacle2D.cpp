/**
 * @file    Obstacle2D.cpp
 * @author  Jon Woolfrey
 * @email   jonathan.woolfrey@gmail.com
 * @date    September 2025
 * @version 1.0
 * @brief   Represents a 2D obstacle with a pose and a shape.
 *
 * @details An Obstacle2D contains a Pose2D (position + orientation) and a pointer
 *          to a Shape2D object. The Shape2D is defined in its local frame.
 *
 * @copyright (c) 2025 Jon Woolfrey
 *
 * @license   OSCL - Free for non-commercial open-source use only.
 *            Commercial use requires a license.
 */
 
#include <Model/Obstacle2D.h>

namespace RobotLibrary { namespace Model {

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                        Constructor                                             //
////////////////////////////////////////////////////////////////////////////////////////////////////
Obstacle2D::Obstacle2D(std::unique_ptr<Shape2D> shape,
                       const double &mass,
                       const double &inertia,
                       const std::string &name)
: RigidBody2D(mass, inertia, name),
  _shape(std::move(shape))
{
    // Worker bees can leave.
    // Even drones can fly away.
    // The Queen is their slave.
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                        Constructor                                             //
////////////////////////////////////////////////////////////////////////////////////////////////////
double
Obstacle2D::distance_to_surface(const Eigen::Vector2d &point) const
{
    Eigen::Vector2d localPoint = _pose.inverse() * point;                                           // Transform to local frame
    
    return _shape->distance_to_surface(localPoint);                                                 // Get distance to shape surface
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                        Constructor                                             //
////////////////////////////////////////////////////////////////////////////////////////////////////
double
Obstacle2D::distance_to_surface_squared(const Eigen::Vector2d &point) const
{
    Eigen::Vector2d localPoint = _pose.inverse() * point;                                           // Convert to local frame
    
    return _shape->distance_to_surface_squared(localPoint);
}

} } // namespace RobotLibrary::Model
