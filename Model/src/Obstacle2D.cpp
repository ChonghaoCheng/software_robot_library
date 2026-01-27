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
    
    if (_shape == nullptr)
    {
        throw std::invalid_argument("[ERROR] [OBSTACLE 2D] Constructor: "
                                    "Shape argument is a null pointer.\n");
    }
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                    Get geometric information about a point relative to this obstacle           //
////////////////////////////////////////////////////////////////////////////////////////////////////
RobotLibrary::Math::ShapeQuery<2>
Obstacle2D::query_point(const Eigen::Vector2d &referencePoint) const
{
    // Transform point to local frame
    Eigen::Vector2d localPoint = _pose.inverse() * referencePoint;

    // Query shape in local frame
    auto query = _shape->query_point(localPoint);

    // Transform results back to world frame
    query.translationVector = _pose.rotation() * query.translationVector;
    query.unitVector        = _pose.rotation() * query.unitVector;
    query.pointOnSurface    = _pose * query.pointOnSurface;

    return query;
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                  Get a point on the surface                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Vector2d
Obstacle2D::point_on_surface(const Eigen::Vector2d &referencePoint) const
{
    Eigen::Vector2d transformedPoint = _pose.inverse() * referencePoint;                            // Transform point to local frame of shape
    
    return _pose * _shape->point_on_surface(transformedPoint);                                      // Get point on surface and map back to world frame
}

} } // namespace RobotLibrary::Model
