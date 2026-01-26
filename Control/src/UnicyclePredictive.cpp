/**
 * @file    UnicyclePredictive.cpp
 * @author  Jon Woolfrey
 * @email   jonathan.woolfrey@gmail.com
 * @date    May 2025
 * @version 1.0
 * @brief   Source files for the MPC class.
 *
 * @details This class solves the predictive control problem for UNICYCLE robot over a
 *          finite number of steps. It uses quadratic functions for both the final step, and all
 *          intermediate steps. The weighting on the final & intermediate pose errors are used as
 *          constructor arguments, whereas the weighting on intermediate control values is based on
 *          the robot's mass & inertia in the RobotLibary::Model::Unicycle class.
 * 
 * @copyright Copyright (c) 2025 Jon Woolfrey
 * 
 * @license GNU General Public License V3
 * 
 * @see https://github.com/Woolfrey/software_robot_library for more information.
 * @see http://github.com/Woolfrey/software_simple_qp for info on the QP solver.
 */
 
#include <Control/UnicyclePredictive.h>

namespace RobotLibrary { namespace Control {

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                            Constructor                                         //
////////////////////////////////////////////////////////////////////////////////////////////////////
UnicyclePredictive::UnicyclePredictive(RobotLibrary::Model::UnicycleParameters &modelParameters,
                                       RobotLibrary::Control::UnicyclePredictiveParameters &controlParameters)
: UnicycleBase(controlParameters.controlFrequency,
               modelParameters),
 _numberOfRecursions(controlParameters.numberOfRecursions),
 _obstaclePotentialScalar(controlParameters.obstaclePotentialScalar),
 _predictionSteps(controlParameters.predictionSteps),
 _threshold(controlParameters.maximumControlStepNorm)
{    
    // Ensure weighting matrices are positive definite:
    std::string message;
    if (not RobotLibrary::Math::is_positive_definite(controlParameters.poseErrorWeight, message))
    {
        throw std::invalid_argument("[ERROR] [UNICYCLE PREDICTIVE CONTROL] Constructor: "
                                    "Initial pose error weight matrix is not positive definite: " + message);
    }
    
    // Check that the exponent is positive
    if (controlParameters.exponent <= 0.0)
    {
        throw std::invalid_argument("[ERROR] [UNICYCLE PREDICTIVE CONTROL] Constructor: "
                                    "Exponent must be positive but received " + std::to_string(controlParameters.exponent) + ".");
    }
 
    // Resize vectors based on prediction horizon
    int N = _predictionSteps;
    _predictedStates.resize(N+1);
    _poseErrorWeight.resize(N+1);
    
    // Generate gain matrices so that M[N-1] == M, and K[N-1] = K
    double a = controlParameters.exponent;
        
    for (int i = 0; i < N; ++i)
    {
        double s = (1.0 / N) + (1.0 - (1.0 / N)) * (1.0 - std::exp(-a * (i / (N - 1.0)))) / (1.0 - std::exp(-a));
        
        // double s = (1.0 / N) + (1.0 - (1.0 / N)) * (std::exp((a * i)/(N - 1.0)) - 1.0) / (std::exp(a) - 1.0);
        
        _poseErrorWeight[i] = s * controlParameters.poseErrorWeight;
    }
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                        Update the state                                        //
////////////////////////////////////////////////////////////////////////////////////////////////////
void
UnicyclePredictive::update_state(const RobotLibrary::Model::Pose2D &pose,
                                 const Eigen::Vector2d &velocity,
                                 const Eigen::Matrix3d &covariance)
{
    RobotLibrary::Model::Unicycle::update_state(pose, velocity, covariance);                        // Update the underlying model
    
    // Transfer these from the base class so we can access them via index in the
    // backward / forward recursions
    _predictedStates[0].pose       = _pose;
    _predictedStates[0].velocity   = this->velocity();                                              // Need "this->" to refer to the method in the base class
    _predictedStates[0].covariance = _covariance;
    
    // Shift the predicted states backward
    for (int i = 1; i < _predictionSteps - 1; ++i)
    {
        _predictedStates[i] = _predictedStates[i+1];
    }
} 

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                               Solve the trajectory tracking problem                            //
////////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Vector2d
UnicyclePredictive::track_trajectory(const std::vector<RobotLibrary::Model::UnicycleState>  &desiredStates,
                                     const std::vector<std::vector<RobotLibrary::Model::Obstacle2D>> &obstacles)
{
    using namespace Eigen;
    using namespace RobotLibrary::Model;
    
    // Ensure inputs are sound
    if (desiredStates.size() != _predictionSteps + 1)
    {
        throw std::invalid_argument("[ERROR] [UNICYCLE PREDICTIVE CONTROL] track_trajectory(): "
                                    "This controller requires N + 1 = " + std::to_string(_predictionSteps+1) + " "
                                    "desired states for the trajectory tracking, but received "
                                    + std::to_string(desiredStates.size()) + ".");
    }
    else if (obstacles.size() != _predictionSteps+1)
    {
        throw std::invalid_argument("[ERROR] [UNICYCLE PREDICTIVE CONTROL] track_trajectory(): "
                                    "This controller has N + 1 = " + std::to_string(_predictionSteps+1) + " "
                                    "prediction steps, but the obstacle array had "
                                    + std::to_string(obstacles.size()) + " elements.");
    }
    
    // Run the optimisation
    for (int i = 0; i < _numberOfRecursions; ++i)
    { 
        double largestStepChange = 0.0;                                                             // Store largest step change in control for this recursion

        double potentialDivisor = 1.0 * i + 1.0;                                                    // Shrinks potential function with each iteration
        
        Vector3d lagrangeMultipliers;                                                               // This equivalent to a wrench for SE(2)
        
        // Backwards recursions
        for (int j = _predictionSteps; j >= 0; --j)
        {
            if (j == _predictionSteps)                                                              // i.e final configuration
            {
                Pose2D currentPose = _predictedStates[j].pose;                                      // This just makes code shorter
                
                Vector3d potentialGradient = -_poseErrorWeight[j] * currentPose.error(desiredStates[j].pose); // NOTE: Force is -K * e = dP/dx
                
                // Compute force from all obstacles
                for (int k = 0; k < obstacles[j].size(); ++k)
                {
                    Vector2d currentPosition = currentPose.translation();                           // For brevity
                    
                    Vector2d pointOnSurface = obstacles[j][k].point_on_surface(currentPosition);    // This is not necessarily the closest point
                    
                    Vector2d translation = currentPosition - pointOnSurface;                        // Translation FROM the surface TO the robot
                    
                    double distance = translation.norm() - _minimumSafeDistance;                    // Subtract safe distance
   
                    if (distance < 0.0)
                    {
                        throw std::runtime_error("[ERROR] [UNICYCLE PREDICTIVE CONTROL] track_trajectory(): "
                                                 "Collision detected on prediction step " + std::to_string(j+1) + " "
                                                 "with obstacle " + std::to_string(k+1) + ".");
                    }

                    /**************************** Harmonic ***************************************/
                    potentialGradient.head(2) -= (_obstaclePotentialScalar / potentialDivisor)
                                               * translation / (pow(distance,2.0) + 1e-08);
                    /******************************************************************************/
                    
                    /****************************** NON Harmonic **********************************
                    potentialGradient.head(2) -= (_obstaclePotentialScalar / potentialDivisor)
                                               * translation / (pow(distance,3.0) + 1e-06);
                    /*******************************************************************************/
                }

                lagrangeMultipliers = -potentialGradient;                                           // This is a force vector              
            }
            else
            {
                // Variables used in this scope
                double angle               = _predictedStates[j].pose.angle();                      // Makes code shorter
                Pose2D currentPose         = _predictedStates[j].pose;
                Pose2D nextPose            = _predictedStates[j+1].pose;
                Vector2d currentVelocity   = _predictedStates[j].velocity;
                Matrix3d potentialHessian  = _poseErrorWeight[j];                                   // i.e. stiffness matrix K
                Vector3d potentialGradient = -potentialHessian * nextPose.error(desiredStates[j+1].pose); // Error at NEXT step, K * e[j+1]
                
                // Add up effects from obstacles
                _distanceToObstacle.resize(obstacles[j+1].size());                                  // We need this to enforce obstacle constraints
                _unitVector.resize(obstacles[j+1].size());                                          // We need this to enforce obstacle constraints
                
                for (int k = 0; k < obstacles[j+1].size(); ++k)
                {
                    Vector2d robotPosition     = nextPose.translation();
                    Vector2d pointOnSurface    = obstacles[j+1][k].point_on_surface(robotPosition);
                    Vector2d translationVector = robotPosition - pointOnSurface;
                    _distanceToObstacle[k]     = translationVector.norm() - _minimumSafeDistance;
                    _unitVector[k]             = translationVector.normalized();

                    if (_distanceToObstacle[k] < 0.0)
                    {           
                        throw std::runtime_error("[ERROR] [UNICYCLE PREDICTIVE CONTROL] track_trajectory(): "
                                                 "Collision detected with '" + obstacles[j+1][k].name() + "' obstacle "
                                                 " on prediction step " + std::to_string(j+1) + ".");
                    }
                    
                    /******************************** Harmonic ************************************/
                    potentialGradient.head(2) -= (_obstaclePotentialScalar / potentialDivisor)
                                               * translationVector / (pow(_distanceToObstacle[k],2.0) + 1e-08);

                    potentialHessian.block(0,0,2,2) += (_obstaclePotentialScalar / (potentialDivisor * (pow(_distanceToObstacle[k],2.0) + 1e-08)))
                                                     * ( 2.0 * (translationVector * translationVector.transpose()) / (pow(_distanceToObstacle[k],2.0) + 1e-08)   - Matrix2d::Identity());
                    /*******************************************************************************/

                    /******************************** NON Harmonic ********************************
                    potentialGradient.head(2) -= (_obstaclePotentialScalar / potentialDivisor)
                                               * translationVector / (pow(_distanceToObstacle[k],3.0) + 1e-06);

                    potentialHessian.block(0,0,2,2) += (_obstaclePotentialScalar / (potentialDivisor * (pow(_distanceToObstacle[k],3.0) + 1e-06)))
                                                     * ( 3.0 * (translationVector * translationVector.transpose()) / (pow(_distanceToObstacle[k],2.0) + 1e-06)   - Matrix2d::Identity());
                    /*******************************************************************************/
                }
                
                // Compute Newton step and update control input u
                
                Vector3d temp = potentialGradient - lagrangeMultipliers;                            // Need this in a couple of places below
                
                Matrix<double,3,3> dfdx = configuration_jacobian(currentPose, currentVelocity, _controlFrequency); // Partial derivative of kinematics w.r.t configuration x
                
                Matrix<double,3,2> dfdu = control_jacobian(currentPose, _controlFrequency);         // Partial derivative of kinematics w.r.t. control input u
                dfdu(2,1) = 1.0;                                                                    // NOTE: This works better for some reason???
                
                Vector<double,2> dLdu = -_inertiaMatrix * (desiredStates[j].velocity - _predictedStates[j].velocity) + dfdu.transpose() * temp; // Partial derivative of Lagrangian w.r.t. configuration x
              
                Matrix<double,2,3> d2Ldudx = dfdu.transpose() * potentialHessian * dfdx;            // Mixed partial derivatives of Lagrangian w.r.t. control u, configuration x
                d2Ldudx(0,2) += (temp[0] * sin(angle) - temp[1] * cos(angle)) / _controlFrequency;  // This is d^2f/dudx^T * (dp/dx - lambda[i+1])
                                           
                Matrix<double,2,2> d2Ldu2 = _inertiaMatrix + dfdu.transpose() * potentialHessian * dfdu; // Second derivative of Lagrangian w.r.t. control u
                
                // NOTE: This somehow prevents obstacle collision,
                // but it can't be too small, or too large
                d2Ldu2(0,0) += 1e-06;
                d2Ldu2(1,1) += 1e-06;
                
                Vector3d dx = currentPose.error(desiredStates[j].pose);
                
                Vector2d du = -d2Ldu2.llt().solve(dLdu + d2Ldudx * dx);                             // Newton step
                
                // Compute instantaneous actuator limits
                Limits linear, angular;
                compute_control_limits(linear, angular, currentVelocity);            
                
                // Clamp
                du[0] = std::clamp(du[0],  linear.lower - currentVelocity[0],  linear.upper - currentVelocity[0]);
                du[1] = std::clamp(du[1], angular.lower - currentVelocity[1], angular.upper - currentVelocity[1]);
                
                double norm = du.norm();
                
                if (norm > largestStepChange) largestStepChange = norm;
                
                _predictedStates[j].velocity += du;                                                 // Update velocity

                lagrangeMultipliers = -potentialGradient + dfdx.transpose() * lagrangeMultipliers;  // Accrue forces
            }
        }
        
        // Forward recursions
        for (int j = 0; j < _predictionSteps; ++j)
        {       
            _predictedStates[j+1].pose =
            RobotLibrary::Model::Unicycle::predicted_pose(_predictedStates[j].pose,
                                                                   _predictedStates[j].velocity,
                                                                   _controlFrequency);
                                              
            _predictedStates[j+1].covariance =
            RobotLibrary::Model::Unicycle::predicted_covariance(_predictedStates[j].pose,
                                                                         _predictedStates[j].velocity,
                                                                         _predictedStates[j].covariance,
                                                                         _controlFrequency);
        }
        
        if (largestStepChange < _threshold) break;
    }
    
    return _predictedStates[0].velocity;                                                            // Only return the 1sts
}

} } // namespace
