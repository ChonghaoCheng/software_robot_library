/**
 * @file    DifferentialDrivePredictive.cpp
 * @author  Jon Woolfrey
 * @email   jonathan.woolfrey@gmail.com
 * @date    May 2025
 * @version 1.0
 * @brief   Source files for the MPC class.
 *
 * @details This class solves the predictive control problem for differential drive robot over a
 *          finite number of steps. It uses quadratic functions for both the final step, and all
 *          intermediate steps. The weighting on the final & intermediate pose errors are used as
 *          constructor arguments, whereas the weighting on intermediate control values is based on
 *          the robot's mass & inertia in the RobotLibary::Model::DifferentialDrive class.
 * 
 * @copyright Copyright (c) 2025 Jon Woolfrey
 * 
 * @license GNU General Public License V3
 * 
 * @see https://github.com/Woolfrey/software_robot_library for more information.
 * @see http://github.com/Woolfrey/software_simple_qp for info on the QP solver.
 */
 
#include <Control/DifferentialDrivePredictive.h>

namespace RobotLibrary { namespace Control {

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                            Constructor                                         //
////////////////////////////////////////////////////////////////////////////////////////////////////
DifferentialDrivePredictive::DifferentialDrivePredictive(RobotLibrary::Model::DifferentialDriveParameters &modelParameters,
                                                         RobotLibrary::Control::DifferentialDrivePredictiveParameters &controlParameters,
                                                         SolverOptions<double> &solverOptions)
: DifferentialDriveBase(controlParameters.controlFrequency,
                        modelParameters,
                        solverOptions),
 _numberOfRecursions(controlParameters.numberOfRecursions),
 _obstaclePotentialScalar(controlParameters.obstaclePotentialScalar),
 _predictionSteps(controlParameters.predictionSteps),
 _threshold(controlParameters.maximumControlStepNorm)
{
    // Ensure weighting matrices are positive definite:
    std::string message;
    if (not RobotLibrary::Math::is_positive_definite(controlParameters.poseErrorWeight, message))
    {
        throw std::invalid_argument("[ERROR] [DIFFERENTIAL DRIVE MPC] Constructor: "
                                    "Initial pose error weight matrix is not positive definite: " + message);
    }
    
    // Set size of vectors
    _predictedStates.resize(_predictionSteps + 1);                                                  // 1 for current state + N for predicted states
    _poseErrorWeight.resize(_predictionSteps);
    _controlWeight.resize(_predictionSteps);

    // Pose Error Weights (Normalized Exponential)
    double exponent = controlParameters.exponent;

    // Precompute denominator for normalization
    double denominator = 0.0;
    for (int j = 0; j < _predictionSteps; ++j)
    {
        denominator += std::exp(exponent * j);
    }

    // Assign normalized exponential pose weights
    for (int j = 0; j < _predictionSteps; ++j)
    {
        double scalar = std::exp(exponent * j) / denominator;
        _poseErrorWeight[j] = scalar * controlParameters.poseErrorWeight;
    }

    std::vector<double> expWeights(_predictionSteps);

    // Compute exponential profile
    for (int i = 0; i < _predictionSteps; ++i)
    {
        expWeights[i] = std::exp(exponent * i);
    }

    // Anchor at R_0 = M or R_{N-1} = M depending on direction of growth
    double scale = (exponent >= 0.0)
                 ? 1.0 / expWeights.back()   // Ensure R_{N-1} = M
                 : 1.0 / expWeights.front(); // Ensure R_0 = M

    // Assign scaled control effort weights
    for (int i = 0; i < _predictionSteps; ++i)
    {
        double weight = expWeights[i] * scale;
        _controlWeight[i] = weight * _inertiaMatrix;
    }                                                 
}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
 //                                        Update the state                                        //
////////////////////////////////////////////////////////////////////////////////////////////////////
void
DifferentialDrivePredictive::update_state(const RobotLibrary::Model::Pose2D &pose,
                                          const Eigen::Vector2d &velocity,
                                          const Eigen::Matrix3d &covariance)
{
    RobotLibrary::Model::DifferentialDrive::update_state(pose, velocity, covariance);               // Update the underlying model
    
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
DifferentialDrivePredictive::track_trajectory(const std::vector<RobotLibrary::Model::DifferentialDriveState>  &desiredStates,
                                              const std::vector<std::vector<RobotLibrary::Model::Obstacle2D>> &obstacles)
{
    using namespace Eigen;
    using namespace RobotLibrary::Model;
    
    // Ensure inputs are sound
    if (desiredStates.size() != _predictionSteps + 1)
    {
        throw std::invalid_argument("[ERROR] [DIFFERENTIAL DRIVE PREDICTIVE] track_trajectory(): "
                                    "This controller requires N + 1 = " + std::to_string(_predictionSteps+1) + " "
                                    "desired states for the trajectory tracking, but received "
                                    + std::to_string(desiredStates.size()) + ".");
    }
    else if (obstacles.size() != _predictionSteps+1)
    {
        throw std::invalid_argument("[ERROR] [DIFFERENTIAL DRIVE PREDICTIVE] track_trajectory(): "
                                    "This controller has N + 1 = " + std::to_string(_predictionSteps+1) + " "
                                    "prediction steps, but the obstacle array had "
                                    + std::to_string(obstacles.size()) + " elements.");
    }
    
    // Run the optimisations
    for (int i = 0; i < _numberOfRecursions; ++i)
    {   
        double largestStepChange = 0.0;                                                             // Store largest step change in control for this recursion

        double potentialDivisor = 1.0 * i + 1;                                                      // Shrinks potential function with each iteration
        
        Vector3d lagrangeMultipliers;                                                               // This equivalent to a wrench for SE(2)
        
        // Backwards recursions
        for (int j = _predictionSteps; j >= 0; --j)
        {
            if (j == _predictionSteps)
            {
                Pose2D currentPose = _predictedStates[j].pose;
                
                Vector3d potentialGradient = - _poseErrorWeight[j] * currentPose.error(desiredStates[j].pose); // NOTE: Force is K * e = - dP/dx
                
                for (int k = 0; k < obstacles[j].size(); ++k)
                {
                    Vector2d currentPosition = currentPose.translation();
                    
                    Vector2d nearestPoint = obstacles[j][k].point_on_surface(currentPosition);
                    
                    Vector2d r = currentPosition - nearestPoint;
                    
                    double distance = r.norm() - _minimumSafeDistance;
                    
                    if (distance <= 0.0)
                    {
                        throw std::runtime_error("[ERROR] [DIFFERENTIAL DRIVE PREDICTIVE] track_trajectory(): "
                                                 "Collision detected on prediction step " + std::to_string(j+1) + " "
                                                 "with obstacle " + std::to_string(k+1) + ".");
                    }
                    
                    potentialGradient.head(2) += - (_obstaclePotentialScalar / potentialDivisor)* r / (distance * distance + 1e-08);
                }

                lagrangeMultipliers = - potentialGradient;                    
            }
            else
            {
                // Variables used in this scope
                double angle = _predictedStates[j].pose.angle();
                
                Pose2D currentPose = _predictedStates[j].pose;
                Pose2D nextPose    = _predictedStates[j+1].pose;
                
                Vector2d currentVelocity = _predictedStates[j].velocity;
                
                Matrix2d M = _controlWeight[j];
                
                Matrix3d potentialHessian = _poseErrorWeight[j];                                    // i.e. stiffness matrix K
                
                Vector3d potentialGradient = - potentialHessian * nextPose.error(desiredStates[j+1].pose); // Error at NEXT step, K * e[j+1]
                
                // Add up effects from obstacles
                for (int k = 0; k < obstacles[j+1].size(); ++k)
                {
                    Vector2d x = nextPose.translation();
                    Vector2d s = obstacles[j+1][k].point_on_surface(x);
                    Vector2d r = x - s;
                    
                    double distance = r.norm() - _minimumSafeDistance;

                    if (distance <= 0.0)
                    {
                        std::cout << "Point on surface: " << s.transpose() << "\n";
                        std::cout << "Robot position:    " << x.transpose() << "\n";
                        std::cout << "Minimum safe distance: " << _minimumSafeDistance << "\n";
                        std::cout << "Distance: " << distance << "\n\n";
                        
                        throw std::runtime_error("[ERROR] [DIFFERENTIAL DRIVE PREDICTIVE] track_trajectory(): "
                                                 "Collision detected on prediction step " + std::to_string(j+1) + " "
                                                 "with obstacle " + std::to_string(k+1) + ".");
                    }
                    
                    double distanceSquared = distance * distance + 1e-08;                           // Add a tiny error to prevent large numbers

                    potentialGradient.head(2) += - (_obstaclePotentialScalar / potentialDivisor) * r / distanceSquared;
                    
                    potentialHessian.block(0,0,2,2) += (_obstaclePotentialScalar / potentialDivisor) * ( 2 * (r * r.transpose()) / distanceSquared - Matrix2d::Identity()) / distanceSquared;
                }
                
                Vector3d temp = potentialGradient - lagrangeMultipliers;
                
                // Partial derivative of kinematics w.r.t configuration x
                Matrix<double,3,3> dfdx = configuration_jacobian(currentPose, currentVelocity, _controlFrequency);
                
                // Partial derivative of kinematics w.r.t. control input u
                Matrix<double,3,2> dfdu = control_jacobian(currentPose, _controlFrequency);
                dfdu(2,1) = 1.0; // NOTE: This works better for some reason???
                
                // Partial derivative of Lagrangian w.r.t. configuration x
                Vector<double,2> dLdu = - M * (desiredStates[j].velocity - _predictedStates[j].velocity) + dfdu.transpose() * temp;
              
                // Mixed partial derivatives of Lagrangian w.r.t. control u, configuration x
                Matrix<double,2,3> d2Ldudx = dfdu.transpose() * potentialHessian * dfdx;
                d2Ldudx(0,2) += (temp[0] * sin(angle) - temp[1] * cos(angle)) / _controlFrequency;  // This is d^2f/dudx^T * (dp/dx - lambda[i+1])
                                           
                // Second derivative of Lagrangian w.r.t. control u
                Matrix<double,2,2> d2Ldu2 = M + dfdu.transpose() * potentialHessian * dfdu;
                
                // Solve for the optimal step size
                Vector3d dx = currentPose.error(desiredStates[j].pose);
                Vector2d du = -d2Ldu2.ldlt().solve(dLdu + d2Ldudx * dx);
                
                double norm = du.norm();
                
                if (norm > largestStepChange) largestStepChange = norm;
                
                _predictedStates[j].velocity += du;
                
                // Constrain
                Limits linear, angular;
                compute_control_limits(linear, angular, currentVelocity);
                
                _predictedStates[j].velocity[0] = std::clamp(_predictedStates[j].velocity[0],  linear.lower,  linear.upper);
                _predictedStates[j].velocity[1] = std::clamp(_predictedStates[j].velocity[1], angular.lower, angular.upper);

                lagrangeMultipliers = -potentialGradient + dfdx.transpose() * lagrangeMultipliers;
            }
        }
        
        // Forward recursions
        for (int j = 0; j < _predictionSteps; ++j)
        {       
            _predictedStates[j+1].pose =
            RobotLibrary::Model::DifferentialDrive::predicted_pose(_predictedStates[j].pose,
                                                                   _predictedStates[j].velocity,
                                                                   _controlFrequency);
                                              
            _predictedStates[j+1].covariance =
            RobotLibrary::Model::DifferentialDrive::predicted_covariance(_predictedStates[j].pose,
                                                                         _predictedStates[j].velocity,
                                                                         _predictedStates[j].covariance,
                                                                         _controlFrequency);
        }
        
        if (largestStepChange < _threshold) break;
    }
    
    return _predictedStates[0].velocity;                                                            // Only return the 1sts
}

} } // namespace
