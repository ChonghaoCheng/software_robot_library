/**
 * @file    RigidBody2D.h
 * @author  Jon Woolfrey
 * @email   jonathan.woolfrey@gmail.com
 * @date    September 2025
 * @version 1.0
 * @brief   A class for describing the kinematics & dynamics of a single rigid body in 2D.
 *
 * @details This class contains the planar kinematic & dynamic properties of a rigid body.
 *          It can be used in SE(2) controllers or 2D physics simulations.
 *
 * @copyright (c) 2025 Jon Woolfrey
 *
 * @license   OSCL - Free for non-commercial open-source use only.
 *            Commercial use requires a license.
 */

#ifndef RIGID_BODY_2D_H
#define RIGID_BODY_2D_H

#include <Model/Pose2D.h>
#include <Eigen/Core>
#include <string>

namespace RobotLibrary { namespace Model {

/**
 * @brief A class representing a solid object in 2D space.
 */
class RigidBody2D
{
    public:
    
        RigidBody2D() = default;

        RigidBody2D(const double &mass,
                    const double &inertia,
                    const std::string &name = "unnamed");

        /**
         * @brief Get the mass (kg) of this object.
         * @return A positive scalar.
         */
        double
        mass() const { return _mass; }
        
        /**
         * @brief Get the inertia (kg*m^2) of this object.
         * @return A positive scalar.
         */
        double
        inertia() const { return _inertia; }
        
        /**
         * @brief Get the combined mass-inertia matrix.
         * @return A 2x2 matrix of the mass and inertia.
         */
        Eigen::Matrix2d
        inertia_matrix() const { return _inertiaMatrix; }
        
        /**
         * @brief Get the current pose of this object.
         * @return An object representing position and orientation in a 2D plane, i.e. SE(2)
         */
        const Pose2D&
        pose() const { return _pose; }
        
        /**
         * @brief Get the current linear and angular velocity of this object.
         * @return A 3D vector containing linear velocity in x & y, and angular velocity.
         */
        const Eigen::Vector3d&
        twist() const { return _twist; }
        
        /**
         * @brief Get the name of this object.
         * @return A string.
         */
        const std::string&
        name() const { return _name; }

        /**
         * @brief Update the state of this object.
         */
        void
        update_state(const Pose2D &pose,
                     const Eigen::Vector3d &twist);

    protected:
    
        double _mass = 0.0;                                                                         ///< Mass of the object (kg)
       
        double _inertia = 0.0;                                                                      ///< Moment of inertia about the z-axis (kg·m²)
       
        Pose2D _pose;                                                                               ///< Pose of the object in global frame
       
        Eigen::Matrix2d _inertiaMatrix;                                                             ///< Combined mass and inertia
        
        Eigen::Vector3d _twist = Eigen::Vector3d::Zero();                                           ///< Linear and angular velocity (v_x, v_y, ω)
        
        std::string _name = "unnamed";                                                              ///< Name / identifier
};

} } // namespace

#endif

