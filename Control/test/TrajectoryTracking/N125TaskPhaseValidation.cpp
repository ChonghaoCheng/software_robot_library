/**
 * @file N125TaskPhaseValidation.cpp
 * @brief Deterministic exact-phase validation using the production RMPCC residual helper.
 */

#include "RmpccPhaseResidual.h"

#include <Math/MathFunctions.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using RobotLibrary::Control::RmpccPhaseAssociation;
using RobotLibrary::Control::RmpccPhaseResiduals;
using RobotLibrary::Control::RmpccStateVector;

constexpr double kPi = 3.14159265358979323846;
constexpr double kRegularization = 1e-8;
constexpr double kObservableTolerance = 1e-12;

struct ErrorPair { double phase = 0.0; double contour = 0.0; };
struct Stats
{
    double bias = 0.0;
    double rms = 0.0;
    double medianAbsolute = 0.0;
    double p95Absolute = 0.0;
    double maximumAbsolute = 0.0;
};
struct MonteCarloSummary
{
    Stats taskPhase;
    Stats metricPhase;
    Stats taskResidual;
    Stats metricResidual;
    double taskCloserFraction = 0.0;
    double rmsDifferenceCiLow = 0.0;
    double rmsDifferenceCiHigh = 0.0;
    double interiorTaskPhaseRms = 0.0;
    double interiorMetricPhaseRms = 0.0;
    int endpointClampedSamples = 0;
    int unobservableSamples = 0;
    int samples = 0;
};

Eigen::Matrix4d canonical_reference(const double s)
{
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T.block<3,1>(0,3) << 0.25 * s,
        0.04 * std::sin(2.0 * kPi * s),
        0.03 * std::sin(kPi * s);
    const Eigen::Vector3d axis = Eigen::Vector3d(0.3, 0.8, 0.5).normalized();
    T.block<3,3>(0,0) =
        (Eigen::AngleAxisd(0.9 * std::sin(kPi * s), axis)
         * Eigen::AngleAxisd(0.5 * s, Eigen::Vector3d::UnitZ())).toRotationMatrix();
    return T;
}

Eigen::Matrix4d attachment(const double d)
{
    Eigen::Matrix4d C = Eigen::Matrix4d::Identity();
    C(2,3) = d;
    return C;
}

Eigen::Matrix4d reference(const double s, const double d)
{
    return canonical_reference(s) * attachment(d);
}

Eigen::Vector<double,6> tangent(const double s, const double d)
{
    constexpr double h = 1e-6;
    const double lo = std::max(0.0, s - h);
    const double hi = std::min(1.0, s + h);
    return RobotLibrary::Math::se3_logarithm(
        RobotLibrary::Math::se3_inverse(reference(lo, d)) * reference(hi, d))
        / (hi - lo);
}

Eigen::Matrix<double,6,6> metric()
{
    Eigen::Matrix<double,6,6> M = Eigen::Matrix<double,6,6>::Identity();
    M.diagonal().tail<3>().setConstant(0.20 * 0.20);
    return M;
}

RmpccPhaseResiduals residual(const Eigen::Matrix4d &actual, const double s,
                             const double d,
                             const RmpccPhaseAssociation association)
{
    RmpccStateVector state = RmpccStateVector::Zero();
    state(6) = s;
    state.head<6>() = RobotLibrary::Math::se3_logarithm(
        RobotLibrary::Math::se3_inverse(reference(s, d)) * actual);
    const auto transform = [d](const double progress)
    {
        return reference(progress, d);
    };
    const auto derivative = [d](const double progress)
    {
        return tangent(progress, d);
    };
    // Exercise both production entry points. The linearized value must retain
    // the exact residual at the nominal state.
    const auto linearized = RobotLibrary::Control::rmpcc_linearize_phase_residuals(
        state, metric(), association, kRegularization, kObservableTolerance,
        1e-6, transform, derivative);
    const auto direct = RobotLibrary::Control::rmpcc_phase_residuals(
        state, metric(), association, kRegularization, kObservableTolerance,
        transform, derivative);
    if((linearized.residual.contour - direct.contour).norm() > 1e-14)
        throw std::runtime_error("production residual and linearization disagree");
    return direct;
}

double position_distance_squared(const Eigen::Matrix4d &actual,
                                 const double s, const double d)
{
    return (actual.block<3,1>(0,3) - reference(s, d).block<3,1>(0,3)).squaredNorm();
}

double exact_task_phase(const Eigen::Matrix4d &actual, const double nominal,
                        const double d)
{
    constexpr int grid = 1001;
    constexpr double lower = 0.05;
    constexpr double upper = 0.95;
    int best = 0;
    double bestValue = std::numeric_limits<double>::infinity();
    for(int i = 0; i < grid; ++i)
    {
        const double s = lower + (upper - lower) * i / (grid - 1.0);
        const double value = position_distance_squared(actual, s, d);
        if(value < bestValue - 1e-18
           || (std::abs(value - bestValue) <= 1e-18
               && std::abs(s - nominal)
                    < std::abs(lower + (upper - lower) * best / (grid - 1.0)
                               - nominal)))
        {
            best = i;
            bestValue = value;
        }
    }
    double a = lower + (upper - lower) * std::max(0, best - 1) / (grid - 1.0);
    double b = lower + (upper - lower) * std::min(grid - 1, best + 1) / (grid - 1.0);
    constexpr double phi = 0.6180339887498948482;
    double x1 = b - phi * (b - a);
    double x2 = a + phi * (b - a);
    double f1 = position_distance_squared(actual, x1, d);
    double f2 = position_distance_squared(actual, x2, d);
    for(int iteration = 0; iteration < 80; ++iteration)
    {
        if(f1 <= f2)
        {
            b = x2; x2 = x1; f2 = f1;
            x1 = b - phi * (b - a);
            f1 = position_distance_squared(actual, x1, d);
        }
        else
        {
            a = x1; x1 = x2; f1 = f2;
            x2 = a + phi * (b - a);
            f2 = position_distance_squared(actual, x2, d);
        }
    }
    return std::clamp(0.5 * (a + b), lower, upper);
}

double fitted_order(const std::vector<double> &scales,
                    const std::vector<double> &errors)
{
    const double floor = 1e-16;
    double mx = 0.0, my = 0.0;
    for(std::size_t i = 0; i < scales.size(); ++i)
    {
        mx += std::log(scales[i]);
        my += std::log(std::max(errors[i], floor));
    }
    mx /= scales.size(); my /= scales.size();
    double numerator = 0.0, denominator = 0.0;
    for(std::size_t i = 0; i < scales.size(); ++i)
    {
        const double x = std::log(scales[i]) - mx;
        numerator += x * (std::log(std::max(errors[i], floor)) - my);
        denominator += x * x;
    }
    return numerator / denominator;
}

Stats statistics(const std::vector<double> &values)
{
    Stats out;
    if(values.empty()) return out;
    out.bias = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
    for(const double value : values) out.rms += value * value;
    out.rms = std::sqrt(out.rms / values.size());
    std::vector<double> absolute;
    absolute.reserve(values.size());
    for(const double value : values) absolute.push_back(std::abs(value));
    std::sort(absolute.begin(), absolute.end());
    out.medianAbsolute = absolute[absolute.size() / 2];
    out.p95Absolute = absolute[static_cast<std::size_t>(0.95 * (absolute.size() - 1))];
    out.maximumAbsolute = absolute.back();
    return out;
}

Eigen::Vector3d random_ball(std::mt19937_64 &generator, const double radius)
{
    std::normal_distribution<double> normal(0.0, 1.0);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    Eigen::Vector3d direction(normal(generator), normal(generator), normal(generator));
    direction.normalize();
    return radius * std::cbrt(unit(generator)) * direction;
}

void write_stats(std::ostream &stream, const Stats &value)
{
    stream << "{\"bias\":" << value.bias << ",\"rms\":" << value.rms
           << ",\"median_absolute\":" << value.medianAbsolute
           << ",\"p95_absolute\":" << value.p95Absolute
           << ",\"maximum_absolute\":" << value.maximumAbsolute << "}";
}

} // namespace

int main(int argc, char **argv)
{
    try
    {
        if(argc != 2)
            throw std::invalid_argument("usage: n125_task_phase_validation OUTPUT_DIR");
        const std::filesystem::path output(argv[1]);
        std::filesystem::create_directories(output);
        std::ofstream samples(output / "n125_task_phase_exact_validation_samples.csv");
        samples << std::setprecision(17)
                << "offset_m,s,delta_s_exact,delta_s_task,delta_s_metric,"
                   "task_phase_error,metric_phase_error,task_residual_error,"
                   "metric_residual_error,task_closer,endpoint_clamped,observable\n";

        const std::array<double,4> magnitudes{0.001, 0.002, 0.004, 0.008};
        std::vector<double> t1PhaseRms, t1ContourRms;
        for(const double magnitude : magnitudes)
        {
            double phaseSquared = 0.0, contourSquared = 0.0;
            int count = 0;
            for(int point = 0; point < 50; ++point)
            {
                const double s = 0.08 + 0.84 * point / 49.0;
                for(const double sign : {-1.0, 1.0})
                {
                    const double ds = sign * magnitude;
                    const Eigen::Matrix4d actual = reference(s + ds, 0.05);
                    const auto task = residual(actual, s, 0.05,
                        RmpccPhaseAssociation::TaskPointXYZ);
                    phaseSquared += std::pow(task.phaseCorrection - ds, 2);
                    contourSquared += task.contour.squaredNorm();
                    ++count;
                }
            }
            t1PhaseRms.push_back(std::sqrt(phaseSquared / count));
            t1ContourRms.push_back(std::sqrt(contourSquared / count));
        }
        const std::vector<double> scaleVector(magnitudes.begin(), magnitudes.end());
        const double t1PhaseSlope = fitted_order(scaleVector, t1PhaseRms);
        const double t1ContourSlope = fitted_order(scaleVector, t1ContourRms);

        double t2TaskMaximum = 0.0, t2MetricMaximum = 0.0;
        const std::array<double,4> degrees{0.5, 1.0, 2.0, 4.0};
        for(int point = 0; point < 50; ++point)
        {
            const double s = 0.08 + 0.84 * point / 49.0;
            std::array<Eigen::Vector3d,4> axes{
                Eigen::Vector3d::UnitX(), Eigen::Vector3d::UnitY(),
                Eigen::Vector3d::UnitZ(), tangent(s, 0.0).tail<3>().normalized()};
            for(const auto &axis : axes) for(const double degree : degrees)
            {
                Eigen::Vector<double,6> xi = Eigen::Vector<double,6>::Zero();
                xi.tail<3>() = degree * kPi / 180.0 * axis;
                const Eigen::Matrix4d actual = reference(s, 0.0)
                    * RobotLibrary::Math::se3_exponential(xi);
                const auto task = residual(actual, s, 0.0,
                    RmpccPhaseAssociation::TaskPointXYZ);
                const auto metricValue = residual(actual, s, 0.0,
                    RmpccPhaseAssociation::MetricScrew);
                t2TaskMaximum = std::max(t2TaskMaximum, std::abs(task.phaseCorrection));
                t2MetricMaximum = std::max(t2MetricMaximum,
                    std::abs(metricValue.phaseCorrection));
            }
        }

        double leverProjectionError = 0.0;
        std::array<double,3> leverResponse{0.0, 0.0, 0.0};
        std::array<double,3> collinearResponse{0.0, 0.0, 0.0};
        for(std::size_t di = 0; di < 3; ++di)
        {
            const double d = 0.05 * di;
            for(int point = 0; point < 30; ++point)
            {
                const double s = 0.08 + 0.84 * point / 29.0;
                std::array<Eigen::Vector3d,4> axes{
                    Eigen::Vector3d::UnitX(), Eigen::Vector3d::UnitY(),
                    Eigen::Vector3d::UnitZ(), tangent(s, 0.0).tail<3>().normalized()};
                for(std::size_t ai = 0; ai < axes.size(); ++ai)
                {
                    Eigen::Vector<double,6> xi = Eigen::Vector<double,6>::Zero();
                    xi.tail<3>() = 0.01 * axes[ai];
                    const Eigen::Matrix4d actual = canonical_reference(s)
                        * RobotLibrary::Math::se3_exponential(xi) * attachment(d);
                    const auto task = residual(actual, s, d,
                        RmpccPhaseAssociation::TaskPointXYZ);
                    const Eigen::Vector3d dp = actual.block<3,1>(0,3)
                        - reference(s, d).block<3,1>(0,3);
                    const Eigen::Vector3d worldTangent = reference(s, d).block<3,3>(0,0)
                        * tangent(s, d).head<3>();
                    const double expected = worldTangent.dot(dp) / worldTangent.squaredNorm();
                    leverProjectionError = std::max(leverProjectionError,
                        std::abs(task.phaseCorrection - expected));
                    leverResponse[di] = std::max(leverResponse[di],
                        std::abs(task.phaseCorrection));
                    if(ai == 2) collinearResponse[di] = std::max(
                        collinearResponse[di], std::abs(task.phaseCorrection));
                }
            }
        }

        std::array<MonteCarloSummary,3> monteCarlo;
        std::mt19937_64 generator(20260817);
        std::uniform_real_distribution<double> progress(0.05, 0.95);
        for(std::size_t di = 0; di < 3; ++di)
        {
            const double d = 0.05 * di;
            std::vector<double> taskPhase, metricPhase, taskContour, metricContour;
            std::vector<double> interiorTaskPhase, interiorMetricPhase;
            int taskCloser = 0;
            for(int sample = 0; sample < 2000; ++sample)
            {
                const double s = progress(generator);
                Eigen::Vector<double,6> xi;
                xi.head<3>() = random_ball(generator, 0.005);
                xi.tail<3>() = random_ball(generator, 5.0 * kPi / 180.0);
                const Eigen::Matrix4d actual = canonical_reference(s)
                    * RobotLibrary::Math::se3_exponential(xi) * attachment(d);
                const double exactProgress = exact_task_phase(actual, s, d);
                const bool endpointClamped = exactProgress <= 0.0500001
                    || exactProgress >= 0.9499999;
                const double exact = exactProgress - s;
                const auto task = residual(actual, s, d,
                    RmpccPhaseAssociation::TaskPointXYZ);
                const auto metricValue = residual(actual, s, d,
                    RmpccPhaseAssociation::MetricScrew);
                const Eigen::Vector<double,6> exactContour =
                    RobotLibrary::Math::se3_logarithm(
                        RobotLibrary::Math::se3_inverse(reference(s + exact, d)) * actual);
                const double taskError = task.phaseCorrection - exact;
                const double metricError = metricValue.phaseCorrection - exact;
                const double taskContourError = (task.contour - exactContour).norm();
                const double metricContourError = (metricValue.contour - exactContour).norm();
                taskPhase.push_back(taskError); metricPhase.push_back(metricError);
                taskContour.push_back(taskContourError); metricContour.push_back(metricContourError);
                taskCloser += std::abs(taskError) < std::abs(metricError);
                if(endpointClamped) ++monteCarlo[di].endpointClampedSamples;
                else
                {
                    interiorTaskPhase.push_back(taskError);
                    interiorMetricPhase.push_back(metricError);
                }
                if(not task.taskPhaseObservable) ++monteCarlo[di].unobservableSamples;
                samples << d << ',' << s << ',' << exact << ','
                        << task.phaseCorrection << ',' << metricValue.phaseCorrection << ','
                        << taskError << ',' << metricError << ',' << taskContourError << ','
                        << metricContourError << ','
                        << (std::abs(taskError) < std::abs(metricError) ? 1 : 0) << ','
                        << (endpointClamped ? 1 : 0) << ','
                        << (task.taskPhaseObservable ? 1 : 0) << '\n';
            }
            auto &summary = monteCarlo[di];
            summary.taskPhase = statistics(taskPhase);
            summary.metricPhase = statistics(metricPhase);
            summary.taskResidual = statistics(taskContour);
            summary.metricResidual = statistics(metricContour);
            summary.taskCloserFraction = static_cast<double>(taskCloser) / taskPhase.size();
            summary.interiorTaskPhaseRms = statistics(interiorTaskPhase).rms;
            summary.interiorMetricPhaseRms = statistics(interiorMetricPhase).rms;
            summary.samples = static_cast<int>(taskPhase.size());
            std::mt19937_64 bootstrap(20260817 + di);
            std::uniform_int_distribution<std::size_t> select(0, taskPhase.size() - 1);
            std::vector<double> differences;
            for(int replicate = 0; replicate < 1000; ++replicate)
            {
                double taskSum = 0.0, metricSum = 0.0;
                for(std::size_t i = 0; i < taskPhase.size(); ++i)
                {
                    const std::size_t index = select(bootstrap);
                    taskSum += taskPhase[index] * taskPhase[index];
                    metricSum += metricPhase[index] * metricPhase[index];
                }
                differences.push_back(std::sqrt(taskSum / taskPhase.size())
                                      - std::sqrt(metricSum / taskPhase.size()));
            }
            std::sort(differences.begin(), differences.end());
            summary.rmsDifferenceCiLow = differences[24];
            summary.rmsDifferenceCiHigh = differences[974];
        }

        std::vector<double> t5Scales{1.0, 0.5, 0.25, 0.125};
        std::vector<double> t5PhaseRms, t5ContourRms;
        std::mt19937_64 scaleGenerator(12520260817ULL);
        std::vector<std::pair<double,Eigen::Vector<double,6>>> bases;
        for(int i = 0; i < 500; ++i)
        {
            Eigen::Vector<double,6> xi;
            xi.head<3>() = random_ball(scaleGenerator, 0.005);
            xi.tail<3>() = random_ball(scaleGenerator, 5.0 * kPi / 180.0);
            bases.emplace_back(progress(scaleGenerator), xi);
        }
        for(const double epsilon : t5Scales)
        {
            double phaseSum = 0.0, contourSum = 0.0;
            for(const auto &[s, xi] : bases)
            {
                const Eigen::Matrix4d actual = canonical_reference(s)
                    * RobotLibrary::Math::se3_exponential(epsilon * xi) * attachment(0.05);
                const double exact = exact_task_phase(actual, s, 0.05) - s;
                const auto task = residual(actual, s, 0.05,
                    RmpccPhaseAssociation::TaskPointXYZ);
                const Eigen::Vector<double,6> exactContour =
                    RobotLibrary::Math::se3_logarithm(
                        RobotLibrary::Math::se3_inverse(reference(s + exact, 0.05)) * actual);
                phaseSum += std::pow(task.phaseCorrection - exact, 2);
                contourSum += (task.contour - exactContour).squaredNorm();
            }
            t5PhaseRms.push_back(std::sqrt(phaseSum / bases.size()));
            t5ContourRms.push_back(std::sqrt(contourSum / bases.size()));
        }
        const double t5PhaseSlope = fitted_order(t5Scales, t5PhaseRms);
        const double t5ContourSlope = fitted_order(t5Scales, t5ContourRms);
        const bool supported = t1PhaseSlope >= 1.8 && t1ContourSlope >= 1.8
            && t5PhaseSlope >= 1.8 && t5ContourSlope >= 1.8
            && monteCarlo[1].interiorTaskPhaseRms
                < monteCarlo[1].interiorMetricPhaseRms
            && monteCarlo[1].unobservableSamples == 0;

        std::ofstream json(output / "n125_task_phase_exact_validation.json");
        json << std::setprecision(17)
             << "{\n  \"schema\":\"n125_task_phase_exact_validation_v1\",\n"
             << "  \"seed\":20260817,\n"
             << "  \"T1\":{\"task_phase_error_slope\":" << t1PhaseSlope
             << ",\"task_contour_residual_slope\":" << t1ContourSlope << "},\n"
             << "  \"T2\":{\"task_phase_max_abs\":" << t2TaskMaximum
             << ",\"metric_phase_max_abs\":" << t2MetricMaximum << "},\n"
             << "  \"T3\":{\"projection_max_abs_error\":" << leverProjectionError
             << ",\"max_abs_response_by_offset\":[" << leverResponse[0] << ','
             << leverResponse[1] << ',' << leverResponse[2]
             << "],\"collinear_z_response_by_offset\":[" << collinearResponse[0]
             << ',' << collinearResponse[1] << ',' << collinearResponse[2] << "]},\n"
             << "  \"T4\":{\n";
        for(std::size_t i = 0; i < monteCarlo.size(); ++i)
        {
            json << "    \"d" << (i * 50) << "\":{\"samples\":" << monteCarlo[i].samples
                 << ",\"task_phase\":"; write_stats(json, monteCarlo[i].taskPhase);
            json << ",\"metric_phase\":"; write_stats(json, monteCarlo[i].metricPhase);
            json << ",\"task_residual\":"; write_stats(json, monteCarlo[i].taskResidual);
            json << ",\"metric_residual\":"; write_stats(json, monteCarlo[i].metricResidual);
            json << ",\"task_closer_fraction\":" << monteCarlo[i].taskCloserFraction
                 << ",\"endpoint_clamped_samples\":"
                 << monteCarlo[i].endpointClampedSamples
                 << ",\"unobservable_samples\":" << monteCarlo[i].unobservableSamples
                 << ",\"interior_task_phase_rms\":"
                 << monteCarlo[i].interiorTaskPhaseRms
                 << ",\"interior_metric_phase_rms\":"
                 << monteCarlo[i].interiorMetricPhaseRms
                 << ",\"task_minus_metric_rms_bootstrap95\":["
                 << monteCarlo[i].rmsDifferenceCiLow << ','
                 << monteCarlo[i].rmsDifferenceCiHigh << "]}"
                 << (i + 1 == monteCarlo.size() ? "\n" : ",\n");
        }
        json << "  },\n  \"T5\":{\"task_phase_error_order\":" << t5PhaseSlope
             << ",\"task_residual_error_order\":" << t5ContourSlope
             << ",\"scales\":[1,0.5,0.25,0.125],\"phase_rms\":[";
        for(std::size_t i = 0; i < t5PhaseRms.size(); ++i)
            json << (i ? "," : "") << t5PhaseRms[i];
        json << "],\"residual_rms\":[";
        for(std::size_t i = 0; i < t5ContourRms.size(); ++i)
            json << (i ? "," : "") << t5ContourRms[i];
        json << "]},\n  \"preliminary_theoretical_support\":"
             << (supported ? "true" : "false") << "\n}\n";

        std::ofstream yaml(output / "n125_task_phase_exact_validation.yaml");
        yaml << std::setprecision(17)
             << "schema: n125_task_phase_exact_validation_v1\n"
             << "seed: 20260817\n"
             << "pure_phase_task_error_slope: " << t1PhaseSlope << "\n"
             << "pure_phase_contour_residual_slope: " << t1ContourSlope << "\n"
             << "rotation_only_task_phase_max_abs: " << t2TaskMaximum << "\n"
             << "rotation_only_metric_phase_max_abs: " << t2MetricMaximum << "\n"
             << "lever_projection_max_abs_error: " << leverProjectionError << "\n"
             << "scaling_task_phase_error_order: " << t5PhaseSlope << "\n"
             << "scaling_task_residual_error_order: " << t5ContourSlope << "\n"
             << "d050_task_phase_rms: " << monteCarlo[1].taskPhase.rms << "\n"
             << "d050_metric_phase_rms: " << monteCarlo[1].metricPhase.rms << "\n"
             << "preliminary_theoretical_support: " << (supported ? "true" : "false") << "\n";

        std::cout << "T1 phase slope=" << t1PhaseSlope
                  << " contour slope=" << t1ContourSlope
                  << "; T5 phase order=" << t5PhaseSlope
                  << " residual order=" << t5ContourSlope
                  << "; d050 task/metric RMS=" << monteCarlo[1].taskPhase.rms
                  << '/' << monteCarlo[1].metricPhase.rms
                  << "; support=" << supported << '\n';
        return supported ? 0 : 2;
    }
    catch(const std::exception &exception)
    {
        std::cerr << exception.what() << '\n';
        return 1;
    }
}
