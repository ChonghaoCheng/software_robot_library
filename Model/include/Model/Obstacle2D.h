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
                   const std::string &name = "unnamed");

        /**
         * @brief Computes the distance from a point to the obstacle's surface.
         * @param point A 2D vector in the global reference frame.
         * @return The (Euclidean?) distance to the surface of the obstacle.
         */
        double
        distance_to_surface(const Eigen::Vector2d &point) const;

        /**
         * @brief Computes the squared distance from a point to the obstacle's surface.
         * @param point A 2D vector in the global reference frame.
         * @return What you asked for.
         */
        double
        distance_to_surface_squared(const Eigen::Vector2d &point) const;

    private:
    
        std::unique_ptr<Shape2D> _shape;                                                            ///< Shape in local frame
};

}} // namespace RobotLibrary::Model

#endif
