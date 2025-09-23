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
 //                                  Get a point on the surface                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Vector2d
Obstacle2D::point_on_surface(const Eigen::Vector2d &referencePoint) const
{
    Eigen::Vector2d transformedPoint = _pose.inverse() * referencePoint;
    
    return _pose * _shape->point_on_surface(transformedPoint);
}

} } // namespace RobotLibrary::Model
