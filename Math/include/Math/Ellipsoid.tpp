/**
 * @file    Ellipsoid.tpp
 * @author  Jon Woolfrey
 * @email   jonathan.woolfrey@gmail.com
 * @date    September 2025
 * @version 1.0
 * @brief   Source files for the Ellipsoid class.
 *
 * @details An ellipsoid is defined by its centerpoint c, and matrix A such that, for any point p on
 *          it surface it satisfies (p - c)^T A^{-1} (p - c) = 1.
 *
 * @copyright (c) 2025 Jon Woolfrey
 *
 * @license   OSCL - Free for non-commercial open-source use only.
 *            Commercial use requires a license.
 * 
 * @see https://github.com/Woolfrey/software_robot_library for more information.
 */

namespace RobotLibrary { namespace Math {

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                          Constructor using the fundamental shape matrix                        //
////////////////////////////////////////////////////////////////////////////////////////////////////
template <unsigned int Dim>
Ellipsoid<Dim>::Ellipsoid(const Eigen::Matrix<double,Dim,Dim> &shapeMatrix)
{
    _shapeMatrix = shapeMatrix;
    
    _LLT = shapeMatrix.llt();
         
    if (_LLT.info() != Eigen::Success)
    {
        throw std::runtime_error("[ERROR] [ELLIPSOID] Constructor: "
                                 "Shape matrix is not positive definite; Cholesky decomposition failed.");
    }
}
        
  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                 Compute the distance to a point                                //
////////////////////////////////////////////////////////////////////////////////////////////////////
template <unsigned int Dim>
double
Ellipsoid<Dim>::distance_to_surface(const Eigen::Vector<double, Dim> &point) const
{
    return sqrt(abs(distance_to_surface_squared(point)));
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                              Compute the squared distance to a point                           //
////////////////////////////////////////////////////////////////////////////////////////////////////
template <unsigned int Dim>
double
Ellipsoid<Dim>::distance_to_surface_squared(const Eigen::Vector<double, Dim> &point) const
{
    return point.dot(_LLT.solve(point)) - 1.0;
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                  Inverse of the shape matrix multiplied by displacement from center            //
////////////////////////////////////////////////////////////////////////////////////////////////////
template <unsigned int Dim>
Eigen::Vector<double, Dim>
Ellipsoid<Dim>::inverse_shape_transformed_vector(const Eigen::Vector<double,Dim> &point) const
{
    return _LLT.solve(point);
}

} } // namespace
