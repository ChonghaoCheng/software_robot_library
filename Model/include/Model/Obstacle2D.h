/**
 * @file    Obstacle2D.h
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

#ifndef OBSTACLE_2D_H
#define OBSTACLE_2D_H

#include <Math/Shape.h>
#include <Model/Pose2D.h>
#include <Model/RigidBody2D.h>
#include <Model/Obstacle2D.h>

#include <memory>

namespace RobotLibrary { namespace Model {

class Obstacle2D : public RigidBody2D
{
    public:

        using Shape2D = RobotLibrary::Math::Shape<2>;

        /**
         * @brief Constructor
         * @param pose The global pose of the obstacle.
         * @param shape A unique pointer to a Shape2D object.
         * @note We need a unique pointer to enable polymorphism.
         */
        Obstacle2D(std::unique_ptr<Shape2D> shape,
                   const double &mass = 1.0,
                   const double &inertia = 1.0,
                   const std::string &name = "");

        /**
         * @brief Get useful geometric information for collision avoidance.
         * @param referencePoint An external reference point.
         * @return A data structure (see Math/include/Math/DataStructures.h)
         */
        RobotLibrary::Math::ShapeQuery<2>
        query_point(const Eigen::Vector2d &referencePoint) const;
        
        /**
         * @brief Get a point on the surface of the obstacle.
         * @note May not be the actual closest point, depending on the underlying shape.
         * @param referencePoint An external point used in the computation.
         * @return A 2D position vector.
         */
        Eigen::Vector2d
        point_on_surface(const Eigen::Vector2d &referencePoint) const;
        
    private:
    
        std::unique_ptr<Shape2D> _shape;                                                            ///< Shape in local frame
};

}} // namespace RobotLibrary::Model

#endif
