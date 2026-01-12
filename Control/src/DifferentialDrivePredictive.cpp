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
        throw std::invalid_argument("[ERROR] [DIFFERENTIAL DRIVE PREDICTIVE] Constructor: "
                                    "Initial pose error weight matrix is not positive definite: " + message);
    }
    
    // Check that the exponent is positive
    if (controlParameters.exponent <= 0.0)
    {
        throw std::invalid_argument("[ERROR] [DIFFERENTIAL DRIVE PREDICTIVE] Constructor: "
                                    "Exponent must be positive but received " + std::to_string(controlParameters.exponent) + ".");
    }
 
    // Resize vectors based on prediction horizon
    int N = _predictionSteps;
    _predictedStates.resize(N+1);
    _poseErrorWeight.resize(N);
    _controlWeight.resize(N);
    
    // Generate gain matrices so that M[N-1] == M, and K[N-1] = K
    double a = controlParameters.exponent;
        
    for (int i = 0; i <= N; ++i)
    {
        double s = (1.0 / N) + (1.0 - (1.0 / N)) * (1.0 - std::exp(-a * (i / (N - 1.0)))) / (1.0 - std::exp(-a));
        
        // double s = (1.0 / N) + (1.0 - (1.0 / N)) * (std::exp((a * i)/(N - 1.0)) - 1.0) / (std::exp(a) - 1.0);
        
        if (i < N)
        {
            _poseErrorWeight[i] = s * controlParameters.poseErrorWeight;
              _controlWeight[i] = s * _inertiaMatrix;
        }
        else
        {
            _finalPoseErrorWeight = s * controlParameters.poseErrorWeight;
        }
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
    
    // Run the optimisation
    for (int i = 0; i < _numberOfRecursions; ++i)
    {   
        double largestStepChange = 0.0;                                                             // Store largest step change in control for this recursion

        double potentialDivisor = 2.0 * i + 1.0;                                                    // Shrinks potential function with each iteration
        
        Vector3d lagrangeMultipliers;                                                               // This equivalent to a wrench for SE(2)
        
        // Backwards recursions
        for (int j = _predictionSteps; j >= 0; --j)
        {
            if (j == _predictionSteps)
            {
                Pose2D currentPose = _predictedStates[j].pose;                                      // This just makes code shorter
                
                Vector3d potentialGradient = - _poseErrorWeight[j] * currentPose.error(desiredStates[j].pose); // NOTE: Force is K * e = - dP/dx
                
                _distanceToObstacle.resize(obstacles[j].size());                                    // Store distance to every obstacle
                
                // Compute force from all obstacles
                for (int k = 0; k < obstacles[j].size(); ++k)
                {
                    Vector2d currentPosition = currentPose.translation();                           // For brevity
                    
                    Vector2d pointOnSurface = obstacles[j][k].pose().translation() + obstacles[j][k].point_on_surface(currentPosition);    // This is not necessarily the closest point
                    
                    Vector2d translation = currentPosition - pointOnSurface;                        // Translation FROM the surface TO the robot
                    
                    Vector2d robotToCentre = obstacles[j][k].pose().translation() - currentPosition; // Distance from robot to obstacle centre

                    Vector2d pointToCentre = obstacles[j][k].pose().translation() - pointOnSurface; //Distance of Surface point to obstacle surface
                    
                    double vectorDir = robotToCentre.norm() - pointToCentre.norm()>0 ? 1.0 :-1.0;

                    
                    double distance = vectorDir * (translation.norm() - _minimumSafeDistance);                    // Store this so we can use it later
                               
                    if (distance <= 0.0)
                    {

                        std::cout << "Point on surface: " << pointOnSurface.transpose() << "\n";
                        std::cout << "Robot position:    " << currentPose.translation().transpose() << "\n";
                        throw std::runtime_error("[ERROR] [DIFFERENTIAL DRIVE PREDICTIVE] track_trajectory(): "
                                                "Collision detected on prediction step " + std::to_string(j+1) + " "
                                                "with obstacle " + std::to_string(k+1) + ".");
                    }
                    
                    potentialGradient.head(2) -= (_obstaclePotentialScalar / potentialDivisor) * translation / (distance * distance + 1e-08);
                }

                lagrangeMultipliers = -potentialGradient;                
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
                
                _distanceToObstacle.resize(obstacles[j+1].size());
                
                _unitVector.resize(obstacles[j+1].size());
                
                // Add up effects from obstacles

                for (int k = 0; k < obstacles[j+1].size(); ++k)
                {
                    Vector2d robotPosition     = nextPose.translation();
                    Vector2d pointOnSurface    =  obstacles[j+1][k].pose().translation() + obstacles[j+1][k].point_on_surface(robotPosition);
                    Vector2d translationVector = robotPosition - pointOnSurface;
                    Vector2d pointToCentre = obstacles[j+1][k].pose().translation() - pointOnSurface;
                    Vector2d robotToCentre = obstacles[j+1][k].pose().translation() - robotPosition;
                    double vectorDir = robotToCentre.norm() - pointToCentre.norm()>0 ? 1.0 :-1.0;
                    
                    _unitVector[k] = translationVector.normalized();
                    
                    _distanceToObstacle[k] = vectorDir *(translationVector.norm() - _minimumSafeDistance);

                    if (_distanceToObstacle[k] <= 0.0)
                    {                       
                        throw std::runtime_error("[ERROR] [DIFFERENTIAL DRIVE PREDICTIVE] track_trajectory(): "
                                                "Collision detected on prediction step " + std::to_string(j+1) + " "
                                                "with obstacle " + std::to_string(k+1) + ".");
                    }
                    
                    double distanceSquared = _distanceToObstacle[k] * _distanceToObstacle[k] + 1e-08; // Add a tiny error to prevent large numbers

                    potentialGradient.head(2) -= (_obstaclePotentialScalar / potentialDivisor) * translationVector / distanceSquared;
                    
                    potentialHessian.block(0,0,2,2) += (_obstaclePotentialScalar / potentialDivisor) * ( 2 * (translationVector * translationVector.transpose()) / distanceSquared - Matrix2d::Identity()) / distanceSquared;
                }
                
                // Compute Newton and update control input u
                
                Vector3d temp = potentialGradient - lagrangeMultipliers;
                
                Matrix<double,3,3> dfdx = configuration_jacobian(currentPose, currentVelocity, _controlFrequency); // Partial derivative of kinematics w.r.t configuration x
                
                Matrix<double,3,2> dfdu = control_jacobian(currentPose, _controlFrequency);         // Partial derivative of kinematics w.r.t. control input u
                dfdu(2,1) = 1.0;                                                                    // NOTE: This works better for some reason???
                
                Vector<double,2> dLdu = - M * (desiredStates[j].velocity - _predictedStates[j].velocity) + dfdu.transpose() * temp; // Partial derivative of Lagrangian w.r.t. configuration x
              
                Matrix<double,2,3> d2Ldudx = dfdu.transpose() * potentialHessian * dfdx;            // Mixed partial derivatives of Lagrangian w.r.t. control u, configuration x
                d2Ldudx(0,2) += (temp[0] * sin(angle) - temp[1] * cos(angle)) / _controlFrequency;  // This is d^2f/dudx^T * (dp/dx - lambda[i+1])
                                           
                Matrix<double,2,2> d2Ldu2 = M + dfdu.transpose() * potentialHessian * dfdu;         // Second derivative of Lagrangian w.r.t. control u
                
                Vector3d dx = currentPose.error(desiredStates[j].pose);                             // Solve for the optimal step size
                
                Vector2d du = -d2Ldu2.llt().solve(dLdu + d2Ldudx * dx);                             // Newton step
                
                // Scale du so we do not violate any boundaries
                Vector2d blah = (dfdu * du).head(2);                                                // Change in configuration due to change in control
                
                double alpha = 1.0;
                
                for (int k = 0; k < obstacles[j].size(); ++k)
                { 
                    double ratio = _distanceToObstacle[k] / blah.dot(_unitVector[k]);
                    
                    if (ratio > 0.0 and ratio < alpha) alpha = 0.99 * ratio;
                }
                
                double norm = alpha * du.norm();
                
                if (norm > largestStepChange) largestStepChange = norm;
                
                _predictedStates[j].velocity += alpha * du;
                
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


