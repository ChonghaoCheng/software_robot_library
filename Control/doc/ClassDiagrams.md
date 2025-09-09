# Control Class Diagrams

[🔙 Back to Control](../README.md)

[🔙 Back to the Foyer](../../README.md)

#### 🧭 Navigation:
- [Data Structures](#data-structures)
- [SerialLink](#seriallink)
   - [SerialLinkBase](#seriallinkbase)
   - [SerialLinkImpedance](#seriallinkimpedance)
   - [SerialLinkKinematic](#seriallinkkinematic)
 
## Data Structures

```mermaid
classDiagram
    direction TB

    class SerialLinkParameters {
        +maxJointAcceleration: double
        +minManipulability: double
        +controlFrequency: uint
        +cartesianPoseGain: Matrix6x6
        +cartesianVelocityGain: Matrix6x6
        +jointPositionGains: vector<double>
        +jointVelocityGains: vector<double>
        +qpsolver: SolverOptions<double>
    }

    class DifferentialDriveFeedbackParameters {
        +controlFrequency: double
        +minimumSafeDistance: double
        +orientationGain: double
        +xPositionGain: double
        +yPositionGain: double
    }

    class DifferentialDrivePredictiveParameters {
        +controlFrequency: double
        +exponent: double
        +maximumControlStepNorm: double
        +minimumSafeDistance: double
        +numberOfRecursions: uint
        +predictionSteps: uint
        +poseErrorWeight: Matrix3x3
    }

    %% Optional: show “used by” relationships
    SerialLinkParameters --> QPSolver : contains options
```

[🔝 Back to Top](#control-class-diagrams)

## SerialLink

### SerialLinkBase

```mermaid
classDiagram
    direction TB

    class QPSolver~double~ {
    }

    class SerialLinkBase {
        <<abstract>>
        +SerialLinkBase(model: shared_ptr<KinematicTree>, endpointName: string, parameters: SerialLinkParameters)
        +resolve_endpoint_motion(endpointMotion: Vector6d) VectorXd*
        +resolve_endpoint_twist(twist: Vector6d) VectorXd*
        +track_endpoint_trajectory(desiredPose: Pose, desiredVel: Vector6d, desiredAcc: Vector6d) VectorXd*
        +track_joint_trajectory(desiredPos: VectorXd, desiredVel: VectorXd, desiredAcc: VectorXd) VectorXd*
        +manipulability() double
        +manipulability_gradient() VectorXd
        +endpoint_pose() Pose
        +pose_error(desired: Pose) Vector6d
        +position_error() double
        +orientation_error() double
        +endpoint_velocity() Vector6d
        +jacobian() Matrix6xN
        +update() void
        +set_redundant_task(task: VectorXd) bool
        +is_singular() bool
        +model() shared_ptr<KinematicTree>
        +frequency() double
        
        -_redundantTaskSet: bool
        -_controlFrequency: double
        -_manipulability: double
        -_maxJointAcceleration: double
        -_minManipulability: double
        -_orientationError: double
        -_positionError: double
        -_cartesianPoseGain: Matrix6x6
        -_cartesianVelocityGain: Matrix6x6
        -_jacobianMatrix: Matrix6xN
        -_forceEllipsoid: Matrix6x6
        -_constraintMatrix: Matrix
        -_constraintVector: VectorXd
        -_redundantTask: VectorXd
        -_model: shared_ptr<KinematicTree>
        -_jointPositionGains: vector<double>
        -_jointVelocityGains: vector<double>
        -_endpointPose: Pose
        -_endpointFrame: ReferenceFrame*
        
        #compute_control_limits(jointNumber: uint) Limits*
    }

    SerialLinkBase --|> QPSolver
    SerialLinkBase --> KinematicTree : uses
    SerialLinkBase --> Pose : uses
    SerialLinkBase --> ReferenceFrame : uses
```

[🔝 Back to Top](#control-class-diagrams)

### SerialLinkImpedance

```mermaid
classDiagram
    direction TB

    class SerialLinkBase {
        <<abstract>>
    }

    class SerialLinkImpedance {

        +SerialLinkImpedance(model: shared_ptr<KinematicTree>, endpointName: string, parameters: SerialLinkParameters)
        +map_wrench_to_torque(wrench: Vector6d) VectorXd
        +set_desired_configuration(configuration: VectorXd) bool
        +resolve_endpoint_motion(endpointMotion: Vector6d) VectorXd
        +resolve_endpoint_twist(twist: Vector6d) VectorXd
        +track_endpoint_trajectory(desiredPose: Pose, desiredVel: Vector6d, desiredAcc: Vector6d) VectorXd
        +track_joint_trajectory(desiredPos: VectorXd, desiredVel: VectorXd, desiredAcc: VectorXd) VectorXd

        -_desiredConfiguration: VectorXd
    }

    SerialLinkImpedance --|> SerialLinkBase
```

[🔝 Back to Top](#control-class-diagrams)

### SerialLinkKinematic

```mermaid
classDiagram
    direction TB

    class SerialLinkBase {
        <<abstract>>
    }

    class SerialLinkKinematic {

        +SerialLinkKinematic(model: shared_ptr<KinematicTree>, endpointName: string, parameters: SerialLinkParameters)
        +resolve_endpoint_motion(endpointMotion: Vector6d) VectorXd
        +resolve_endpoint_twist(twist: Vector6d) VectorXd
        +track_endpoint_trajectory(desiredPose: Pose, desiredVel: Vector6d, desiredAcc: Vector6d) VectorXd
        +track_joint_trajectory(desiredPos: VectorXd, desiredVel: VectorXd, desiredAcc: VectorXd) VectorXd

        #compute_control_limits(jointNumber: uint) Limits
    }

    SerialLinkKinematic --|> SerialLinkBase
```

[🔝 Back to Top](#control-class-diagrams)

