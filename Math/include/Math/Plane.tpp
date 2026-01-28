/**
 * @file    Plane.tpp
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
 
 
#include <cassert>

namespace RobotLibrary { namespace Math {

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                           Constructor                                          //
////////////////////////////////////////////////////////////////////////////////////////////////////
template <unsigned int AmbientDim, unsigned int SubDim>
Plane<AmbientDim, SubDim>::Plane(const BasisType &basis)
{
    // Check arguments are sound
    if (basis.rows() <= basis.cols())
    {
        throw std::invalid_argument("[ERROR] [PLANE] Constructor: " 
                                    "Basis must have more rows than columns: "
                                    + std::to_string(basis.cols()) + " /> " + basis.rows() + ".");
    }
    
    if (basis.cols() == 1) this->_type = "line";
    else                   this->_type = "plane";
    
    // Orthonormalize using Eigen's ColPivHouseholderQR
    Eigen::ColPivHouseholderQR<BasisType> qr(basis);
    _basis = qr.householderQ().leftCols(SubDim);
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                  Compute a point on the surface                                //
////////////////////////////////////////////////////////////////////////////////////////////////////
template <unsigned int AmbientDim, unsigned int SubDim>
Eigen::Vector<double, AmbientDim>
Plane<AmbientDim, SubDim>::point_on_surface(const Eigen::Vector<double,AmbientDim> &referencePoint) const
{
    return _basis * (_basis.transpose() * referencePoint); 
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                    Compute geometry of a reference point relative to this shape                //
////////////////////////////////////////////////////////////////////////////////////////////////////
template <unsigned int AmbientDim, unsigned int SubDim>
ShapeQuery<AmbientDim>
Plane<AmbientDim, SubDim>::query_point(const Eigen::Matrix<double, AmbientDim, 1> &referencePoint) const
{
    ShapeQuery<AmbientDim> query;                                                                   // We want to return this

    query.pointOnSurface = point_on_surface(referencePoint);                                        // Speaks for itself
    
    query.translationVector = referencePoint - query.pointOnSurface;                                // Fine so far

    // Oriented distance
    if constexpr (AmbientDim == 2 and SubDim == 1)
    {
        Eigen::Vector2d direction = _basis.col(0);                                                  // Line direction

        Eigen::Vector2d perpendicularLine(-direction.y(), direction.x());
        
        query.signedDistance = perpendicularLine.dot(query.translationVector);
    }
    else if constexpr (AmbientDim == 3 and SubDim == 2)
    {
        Eigen::Vector3d perpendicular = _basis.col(0).cross(_basis.col(1));
        perpendicular.normalize();
        query.signedDistance = perpendicular.dot(query.translationVector);
    }
    else
    {
        query.signedDistance = query.translationVector.norm();
    }
    
    double norm = query.translationVector.norm();
    
    query.unitVector = (norm > 1e-12)
                     ? (query.signedDistance >= 0.0 ? query.translationVector / norm
                                                     : -query.translationVector / norm)
                     : Eigen::Matrix<double, AmbientDim, 1>::Zero();

    return query;
}

}} // namespace
