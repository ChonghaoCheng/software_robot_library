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
    _shapeMatrix = shapeMatrix;                                                                     // Save internally
    
    _LLT = shapeMatrix.llt();                                                                       // Pre-compute decomposition to save time
         
    if (not _LLT.info() == Eigen::Success)
    {
        throw std::runtime_error("[ERROR] [ELLIPSOID] Constructor: "
                                 "Shape matrix is not positive definite; Cholesky decomposition failed.");
    }
}
        
  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                       Get a point on the surface along a ray to the center                     //
////////////////////////////////////////////////////////////////////////////////////////////////////
template <unsigned int Dim>
Eigen::Vector<double,Dim>
Ellipsoid<Dim>::point_on_surface(const Eigen::Vector<double,Dim> &referencePoint) const
{
    // Ellipsoid: E = { x | x^T A^{-1} x = 1 } (centered at origin)
    // Given a reference direction r, find the intersection of the ray x = λ r
    // with the ellipsoid surface.
    // Solve (λ r)^T A^{-1} (λ r) = 1  ⇒  λ = 1 / sqrt(r^T A^{-1} r)
    // => p = r / sqrt(r^T A^{-1} r)

    return referencePoint / sqrt(referencePoint.dot(_LLT.solve(referencePoint)));
}

}} // namespace
