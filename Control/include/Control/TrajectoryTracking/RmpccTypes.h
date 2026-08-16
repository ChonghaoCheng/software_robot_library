#ifndef CONTROL_TRAJECTORY_TRACKING_RMPCC_TYPES_H
#define CONTROL_TRAJECTORY_TRACKING_RMPCC_TYPES_H

namespace RobotLibrary { namespace Control {

enum class RmpccPredictorGeometry
{
    ExactSE3,
    AdditiveLieAlgebra
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

enum class RmpccResidualLinearization
{
    FrozenProjector,
    FullResidualJacobian
};

} } // namespace RobotLibrary::Control

#endif
