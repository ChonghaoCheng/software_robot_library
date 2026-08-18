#ifndef CONTROL_TRAJECTORY_TRACKING_RMPCC_TYPES_H
#define CONTROL_TRAJECTORY_TRACKING_RMPCC_TYPES_H

namespace RobotLibrary { namespace Control {

enum class RmpccPredictorGeometry
{
    ExactSE3,
    AdditiveLieAlgebra
};

enum class RmpccReferenceMotion
{
    LegacyTangentProduct,
    StageConsistent
};

enum class RmpccLagGeometry
{
    FullScrew,
    SplitTranslationRotation,
    TranslationOnly,
    RotationOnly
};

enum class RmpccObjectiveGeometry
{
    FullScrewSE3,
    DecoupledCartesianSO3
};

enum class RmpccContourResidualGeometry
{
    LocalUnifiedSE3,
    ExactAssociatedUnifiedSE3,
    ScheduledDecoupledCartesianSO3,
    AssociatedDecoupledCartesianSO3
};

enum class RmpccResidualLinearization
{
    FrozenProjector,
    FullResidualJacobian
};

enum class RmpccPhaseAssociation
{
    MetricScrew,
    TaskPointXYZ,
    TaskPoseFeature
};

enum class RmpccLagPenalty
{
    PhaseInducedPoseVector,
    ScalarTaskDistance,
    ScalarPosePathArc
};

} } // namespace RobotLibrary::Control

#endif
