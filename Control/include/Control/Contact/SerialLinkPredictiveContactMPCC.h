/**
 * @file SerialLinkPredictiveContactMPCC.h
 * @brief Seven-input predictive contact MPCC using the E05 M3 plant.
 */

#ifndef SERIAL_LINK_PREDICTIVE_CONTACT_MPCC_H
#define SERIAL_LINK_PREDICTIVE_CONTACT_MPCC_H

#include <Control/Contact/PredictiveContactModel.h>
#include <Control/TrajectoryTracking/SerialLinkMPCC.h>

namespace RobotLibrary { namespace Control {

enum class PredictiveContactMode
{
    Disabled,
    Shadow,
    Active
};

struct PredictiveContactMpccParameters
{
    PredictiveContactMode mode = PredictiveContactMode::Disabled;
    double identifiedDt = 0.002;
    double stiffness = 21303.75539503847;
    double targetForce = 2.5;
    double forceWeight = 1.0;
    double minimumForce = 1.8;
    double maximumForce = 3.2;
    double maximumPenetrationIncrement = 4.0e-6;
    double forceSlackWeight = 1.0e5;
    double penetrationSlackWeight = 1.0e12;
    double maximumForceSlack = 5.0;
    double maximumPenetrationSlack = 1.0e-3;
    double tangentPositionWeight = 260.0;
    double pathLagPositionWeight = 35.0;
    Eigen::Vector3d compressionDirectionParent{0.0, -1.0, 0.0};
    Eigen::Vector3d contactOffsetEndpoint{
        0.0005265895, 0.0999391081, -0.0399956612};
};

struct PredictiveContactMpccDiagnostics
{
    PredictiveContactMode mode = PredictiveContactMode::Disabled;
    bool forceValid = false;
    bool forceObjectiveActive = false;
    double measuredForce = 0.0;
    double stiffness = 0.0;
    double forceCost = 0.0;
    double forceGradientNorm = 0.0;
    double maximumForceSlack = 0.0;
    double maximumPenetrationSlack = 0.0;
    double optimizedRelativeNormalVelocityStage0 = 0.0;
    Eigen::VectorXd predictedRelativeNormalVelocity;
    Eigen::VectorXd predictedPenetrationIncrement;
    Eigen::VectorXd predictedForce;
    Eigen::MatrixXd forceMap;
    Eigen::VectorXd forceOffset;
};

class SerialLinkPredictiveContactMPCC : public SerialLinkMPCC
{
    public:
        using SerialLinkMPCC::SerialLinkMPCC;

        void set_predictive_contact_parameters(
            const PredictiveContactMpccParameters &parameters);

        /** A stale/invalid cycle must explicitly pass valid=false. */
        void set_force_measurement(double measuredForce, bool valid);

        const PredictiveContactMpccParameters &predictive_contact_parameters() const
        {
            return _parameters;
        }

        const PredictiveContactMpccDiagnostics &predictive_contact_diagnostics() const
        {
            return _diagnostics;
        }

    protected:
        Eigen::Matrix3d position_error_weight(
            int stage,
            const Eigen::Vector<double,ERROR_DIM> &pathTangent,
            const Eigen::Matrix3d &defaultWeight,
            const Eigen::Matrix3d &predictionRotation,
            const Eigen::Matrix4d &predictedParentTransform) const override;

        void extend_qp_problem(
            const MpccQpExtensionContext &context,
            Eigen::MatrixXd &hessian,
            Eigen::VectorXd &gradient,
            Eigen::MatrixXd &constraintMatrix,
            Eigen::VectorXd &constraintVector,
            Eigen::VectorXd &seed) override;

        void shift_extension_warm_start(
            const MpccQpExtensionContext &context,
            const Eigen::VectorXd &optimum,
            Eigen::VectorXd &shiftedWarmStart) override;

        void on_extended_qp_solution(
            const MpccQpExtensionContext &context,
            const Eigen::VectorXd &optimum) override;

        void reset_additional_virtual_state() override;

    private:
        bool prediction_enabled() const;
        PredictiveContactKinematics make_kinematics(
            const MpccQpExtensionContext &context) const;

        PredictiveContactMpccParameters _parameters;
        PredictiveContactMpccDiagnostics _diagnostics;
        double _measuredForce = 0.0;
        bool _forceValid = false;
        PredictiveContactAffineModel _lastModel;
        Eigen::MatrixXd _lastForceMap;
        Eigen::VectorXd _lastForceConstant;
        int _baseControlDimension = 0;
};

} } // namespace RobotLibrary::Control

#endif
