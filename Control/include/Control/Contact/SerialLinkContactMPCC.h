/**
 * @file SerialLinkContactMPCC.h
 * @brief Coupled path-progress and normal-force predictive control.
 */

#ifndef SERIAL_LINK_CONTACT_MPCC_H
#define SERIAL_LINK_CONTACT_MPCC_H

#include <Control/Contact/ContactDataStructures.h>
#include <Control/TrajectoryTracking/SerialLinkMPCC.h>

namespace RobotLibrary { namespace Control {

/**
 * MPCC with a rigid-motion-consistent board-relative force model. Loss mode
 * keeps the original 7N decision. Constraint mode appends N force slacks.
 * Virtual progress never enters the physical tool-position/force map directly.
 */
class SerialLinkContactMPCC : public SerialLinkMPCC
{
    public:
        using SerialLinkMPCC::SerialLinkMPCC;

        void set_contact_parameters(const ContactParameters &parameters);
        void set_contact_state(const ContactState &state);

        const ContactParameters &contact_parameters() const { return _contactParameters; }
        const ContactMpcDiagnostics &contact_diagnostics() const { return _contactDiagnostics; }

    protected:
        Eigen::Matrix3d
        position_error_weight(
            int stage,
            const Eigen::Vector<double,ERROR_DIM> &pathTangent,
            const Eigen::Matrix3d &defaultWeight,
            const Eigen::Matrix3d &predictionRotation,
            const Eigen::Matrix4d &predictedParentTransform) const override;

        void
        extend_qp_problem(const MpccQpExtensionContext &context,
                          Eigen::MatrixXd &hessian,
                          Eigen::VectorXd &gradient,
                          Eigen::MatrixXd &constraintMatrix,
                          Eigen::VectorXd &constraintVector,
                          Eigen::VectorXd &seed) override;

        void
        shift_extension_warm_start(const MpccQpExtensionContext &context,
                                   const Eigen::VectorXd &optimum,
                                   Eigen::VectorXd &shiftedWarmStart) override;

        void
        on_extended_qp_solution(const MpccQpExtensionContext &context,
                                const Eigen::VectorXd &optimum) override;

        void
        on_twist_resolved(const Eigen::Vector<double,6> &commandedBaseTwist,
                          const Eigen::Vector<double,6> &realizedBaseTwist,
                          double dt) override;

    private:
        bool contact_active() const;

        ContactParameters _contactParameters;
        ContactState _contactState;
        ContactMpcDiagnostics _contactDiagnostics;

        Eigen::RowVectorXd _firstForceMap;
        double _firstForceConstant = 0.0;
        double _currentSignedNormalCoordinate = 0.0;
        Eigen::Vector3d _firstStageNormal = Eigen::Vector3d::UnitZ();
        Eigen::Vector3d _firstStageBoardPosition = Eigen::Vector3d::Zero();
        Eigen::Vector3d _currentToolPosition = Eigen::Vector3d::Zero();
};

} } // namespace RobotLibrary::Control

#endif
