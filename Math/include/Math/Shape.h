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

#include <Math/DataStructures.h>

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
        
        /**
         * @brief Query the geometry of a point relative to this shape.
         * @param referencePoint Used for computing properties.
         * @return A data structure containg signed distance, translation vector, etc.
         */
        virtual
        ShapeQuery<Dim>
        query_point(const VectorType &referencePoint) const = 0;
        
        /**
         * @brief Get the type of shape.
         */
        std::string
        type() const { return _type; }
        
    protected:
        
        std::string _type = "unknown";
};

// Convenience aliases
using Shape2D = Shape<2>;
using Shape3D = Shape<3>;

} } // namespace RobotLibrary::Math

#endif
