/**
 * @file SerialLinkContactAdaptiveMPCC.h
 * @brief Contact-aware MPCC with contact-coordinate rate as a QP input.
 */

#ifndef SERIAL_LINK_CONTACT_ADAPTIVE_MPCC_H
#define SERIAL_LINK_CONTACT_ADAPTIVE_MPCC_H

#include <Control/TrajectoryTracking/SerialLinkMPCC.h>

namespace RobotLibrary { namespace Control {

struct ContactAdaptiveMpccParameters
{
    bool adaptationEnabled = true;
    bool contactCostEnabled = true;
    bool optimizeDeltaRate = true;
    bool normalAdmittanceEnabled = false;
    double targetForce = 2.5;              ///< Positive compression [N].
    double deltaGain = 1.2e-4;             ///< m/(N s).
    double forceDeadband = 0.1;            ///< N.
    double deltaMinimum = -5e-4;           ///< m.
    double deltaMaximum = 5e-4;            ///< m.
    double deltaRateMaximum = 5e-4;        ///< m/s.
    double forceRateWeight = 1e4;
    double deltaRateSmoothWeight = 0.0;
    double normalAdmittanceGain = 1.2e-4;  ///< m/(N s).
    Eigen::Vector3d compressionDirectionInParent{0.0, -1.0, 0.0};
};

struct ContactAdaptiveMpccDiagnostics
{
    double measuredForce = 0.0;
    double targetForce = 0.0;
    double forceError = 0.0;
    double preferredDeltaRate = 0.0;
    double optimizedDeltaRate = 0.0;
    double delta = 0.0;
    double nextDelta = 0.0;
    bool deltaSaturated = false;
    bool deltaRateSaturated = false;
    Eigen::Vector<double,6> firstDeltaTangent =
        Eigen::Vector<double,6>::Zero();
    Eigen::VectorXd optimizedDeltaRates;
    Eigen::VectorXd predictedDelta;
};

/** Signed deadband that preserves the residual magnitude beyond the band. */
double contact_force_deadband(double forceError, double deadband);

/** Preferred contact-coordinate rate from the current measured force only. */
double contact_preferred_delta_rate(
    const ContactAdaptiveMpccParameters &parameters,
    double measuredForce);

/** Body-coordinate derivative dT_r/d_delta = [R_d^T n_c^B; 0]. */
Eigen::Vector<double,6> contact_delta_tangent_body(
    const Eigen::Matrix3d &pathRotationInParent,
    const Eigen::Vector3d &compressionDirectionInParent);

/** Same tangent expressed in the MPCC's fixed local prediction frame. */
Eigen::Vector<double,6> contact_delta_tangent_prediction_frame(
    const Eigen::Matrix3d &predictionRotationInBase,
    const Eigen::Matrix3d &parentRotationInBase,
    const Eigen::Vector3d &compressionDirectionInParent);

class SerialLinkContactAdaptiveMPCC : public SerialLinkMPCC
{
    public:
        using SerialLinkMPCC::SerialLinkMPCC;

        void set_contact_adaptive_parameters(
            const ContactAdaptiveMpccParameters &parameters);

        void set_measured_normal_force(double measuredForce);

        void set_contact_coordinate(double delta);

        /** CA-MPCC-A update: integrate measured-force rate before the MPCC step. */
        void update_contact_coordinate_from_force(double dt);

        double contact_coordinate() const { return _delta; }

        const ContactAdaptiveMpccParameters &contact_adaptive_parameters() const
        {
            return _parameters;
        }

        const ContactAdaptiveMpccDiagnostics &contact_adaptive_diagnostics() const
        {
            return _contactDiagnostics;
        }

    protected:
        RobotLibrary::Model::Pose
        reference_pose_at_progress(double progress) override;

        int stage_control_dimension() const override
        {
            return _parameters.optimizeDeltaRate ? 8 : 7;
        }

        Eigen::MatrixXd
        additional_reference_tangents(
            int stage,
            double progress,
            const Eigen::Matrix3d &predictionRotation,
            const Eigen::Matrix4d &predictedParentTransform) const override;

        void
        configure_additional_stage_inputs(int stage,
                                          Eigen::VectorXd &lower,
                                          Eigen::VectorXd &upper,
                                          Eigen::VectorXd &nominal) const override;

        void
        extend_qp_problem(const MpccQpExtensionContext &context,
                          Eigen::MatrixXd &hessian,
                          Eigen::VectorXd &gradient,
                          Eigen::MatrixXd &constraintMatrix,
                          Eigen::VectorXd &constraintVector,
                          Eigen::VectorXd &seed) override;

        void
        on_extended_qp_solution(const MpccQpExtensionContext &context,
                                const Eigen::VectorXd &optimum) override;

        void reset_additional_virtual_state() override;

        Eigen::Vector<double,6>
        postprocess_base_twist(const Eigen::Vector<double,6> &baseTwist,
                               double dt) override;

    private:
        static constexpr int DELTA_RATE_INDEX = 7;

        ContactAdaptiveMpccParameters _parameters;
        ContactAdaptiveMpccDiagnostics _contactDiagnostics;
        double _measuredForce = 0.0;
        double _delta = 0.0;
        double _lastDeltaRate = 0.0;
};

} } // namespace RobotLibrary::Control

#endif
