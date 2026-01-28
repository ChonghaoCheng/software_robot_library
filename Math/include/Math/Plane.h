/**
 * @file    Plane.h
 * @author  Jon Woolfrey
 * @email   jonathan.woolfrey@gmail.com
 * @date    February 2026
 * @version 1.0
 * @brief   A class that represents a plane in n-dimensional space.
 *
 * @copyright (c) 2026 Jon Woolfrey
 *
 * @license   OSCL - Free for non-commercial open-source use only.
 *            Commercial use requires a license.
 * 
 * @see https://github.com/Woolfrey/software_robot_library for more information.
 */

#ifndef PLANE_H
#define PLANE_H

#include <Math/Shape.h>                                                                             // Base class

#include <Eigen/Dense>                                                                              // For vector algebra

namespace RobotLibrary { namespace Math {

template <unsigned int AmbientDim, unsigned int SubDim>
class Plane : public Shape<Dim>
{
    public:
    
        using VectorType = Eigen::Matrix<double, AmbientDim, 1>;
        using BasisType  = Eigen::Matrix<double, AmbientDim, SubDim>;

        /**
         * @brief Constructor
         * @param basis Columns are spanning vectors of the subspace
         */
        Plane(const BasisType &basis);
        
        /**
         * @brief Get the point on the surface nearest to a reference point.
         * @param referencePoint The point external to the plane.
         * @return A point on the surface.
         */
        VectorType
        point_on_surface(const VectorType &referencePoint) const override;
        
        /**
         * @brief
         * @param
         * @return A data structure containing information like signed distance, translation vector, etc.
         */
        ShapeQuery<Dim>
        query_point(const VectorType &referencePoint) const override;

    private:
    
        BasisType _basis;
};

// Convenient aliases
using Line2D = Plane<2,1>
using Line3D = Plane<3,1>
using Plan2D = Plane<3,2>

}} // namespace
