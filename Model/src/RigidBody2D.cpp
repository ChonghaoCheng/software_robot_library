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
  
#include <Model/RigidBody2D.h>

namespace RobotLibrary { namespace Model {

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                         Constructor                                            //
////////////////////////////////////////////////////////////////////////////////////////////////////
RigidBody2D::RigidBody2D(const double &mass,
                         const double &inertia,
                         const std::string &name)
: _mass(mass),
  _inertia(inertia),
  _name(name)
{
    // Check the input parameters are sound
    if (_mass <= 0.0)
    {
        throw std::invalid_argument("[ERROR] [RIGID BODY 2D] Constructor: "
                                    "Mass must be positive but argument was " + std::to_string(_mass) + ".");
    }
    else if (_inertia < 0.0)
    {
        throw std::invalid_argument("[ERROR] [RIGID BODY 2D] Constructor: "
                                    "Inertia must be positive but argument was " + std::to_string(_inertia) + ".");
    }
    
    _inertiaMatrix << _mass,      0.0,
                        0.0, _inertia;
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                            Update the position, orientation, and speed                         //
////////////////////////////////////////////////////////////////////////////////////////////////////    
void
RigidBody2D::update_state(const Pose2D &pose,
                          const Eigen::Vector3d &twist)
{
    _pose  = pose;
    _twist = twist;
}

} } // namespace RobotLibary::Model
