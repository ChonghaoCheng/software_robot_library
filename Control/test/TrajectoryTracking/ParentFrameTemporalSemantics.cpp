#include <Control/TrajectoryTracking/ParentFrameReferenceMotion.h>
#include <Math/MathFunctions.h>

#include <Eigen/Core>

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Motion = RobotLibrary::Control::CausalParentFrameMotion;
using Status = Motion::UpdateStatus;
using RobotLibrary::Control::predicted_parent_frame_state;
using RobotLibrary::Math::se3_exponential;
using RobotLibrary::Math::se3_inverse;
using RobotLibrary::Math::se3_logarithm;

int failures = 0;

void check(const bool condition, const std::string &description)
{
    if(!condition)
    {
        ++failures;
        std::cerr << "FAIL: " << description << '\n';
    }
}

Eigen::Matrix4d advanced(const Eigen::Matrix4d &origin,
                         const Eigen::Vector<double,6> &bodyTwist,
                         const double elapsed)
{
    return origin * se3_exponential(elapsed * bodyTwist);
}

double pose_error(const Eigen::Matrix4d &actual, const Eigen::Matrix4d &expected)
{
    return se3_logarithm(se3_inverse(expected) * actual).norm();
}

const std::vector<Eigen::Vector<double,6>> &twist_cases()
{
    static const std::vector<Eigen::Vector<double,6>> cases = {
        (Eigen::Vector<double,6>() << .08, -.03, .04, 0., 0., 0.).finished(),
        (Eigen::Vector<double,6>() << 0., 0., 0., .12, -.09, .07).finished(),
        (Eigen::Vector<double,6>() << .05, .02, -.01, -.08, .10, .04).finished()};
    return cases;
}

void t0_normal_constant_twist()
{
    const Eigen::Matrix4d origin = se3_exponential(
        (Eigen::Vector<double,6>() << .03, -.02, .04, .1, -.07, .05).finished());
    for(std::size_t index = 0; index < twist_cases().size(); ++index)
    {
        Motion motion;
        const auto &twist = twist_cases()[index];
        check(motion.update(origin, 1.0, 4, 1.0) == Status::FirstTimestampedSample,
              "T0 first sample status " + std::to_string(index));
        const Eigen::Matrix4d next = advanced(origin, twist, .010);
        check(motion.update(next, 1.010, 4, 1.010) == Status::TimestampAdvanced,
              "T0 advanced status " + std::to_string(index));
        check((motion.body_twist() - twist).norm() < 2e-12,
              "T0 recovered body twist " + std::to_string(index));
        check(pose_error(motion.predicted_pose(20, .002), advanced(next, twist, .040)) < 2e-12,
              "T0 predicted pose " + std::to_string(index));
        const auto state = predicted_parent_frame_state(motion, 20, .002);
        check(state.twistInBase.allFinite() && state.twistInBase.norm() > 0.0,
              "T0 predicted frame twist " + std::to_string(index));
    }
}

void t1_low_rate_source()
{
    Motion motion;
    const auto twist = twist_cases().back();
    Eigen::Matrix4d sample = Eigen::Matrix4d::Identity();
    std::uint64_t advancedCount = 0;
    std::uint64_t noNewCount = 0;
    for(int source = 0; source < 5; ++source)
    {
        const double measurementTime = 2.0 + .010 * source;
        sample = advanced(Eigen::Matrix4d::Identity(), twist, .010 * source);
        for(int evaluation = 0; evaluation < 5; ++evaluation)
        {
            const auto status = motion.update(
                sample, measurementTime, 9, measurementTime + .002 * evaluation);
            advancedCount += status == Status::TimestampAdvanced;
            noNewCount += status == Status::NoNewMeasurement;
            if(source > 0)
            {
                check((motion.body_twist() - twist).norm() < 2e-12,
                      "T1 duplicate does not reset body twist");
                check(motion.prediction_velocity_active(),
                      "T1 prediction remains fresh");
            }
        }
    }
    check(advancedCount == 4, "T1 one finite difference per new source timestamp");
    check(noNewCount == 20, "T1 four no-new statuses per source interval");
    check(motion.duplicate_timestamp_count() == 20, "T1 duplicate counter");
    check(motion.body_twist().allFinite(), "T1 no division by zero");
}

void t2_duplicate_pose_mismatch()
{
    Motion motion;
    const Eigen::Matrix4d accepted = Eigen::Matrix4d::Identity();
    motion.update(accepted, 3.0, 1, 3.0);
    const Eigen::Matrix4d changed = advanced(accepted, twist_cases()[0], .010);
    const auto status = motion.update(changed, 3.0, 1, 3.002);
    check(status == Status::DuplicateTimestampPoseMismatchRejected, "T2 status");
    check(motion.duplicate_pose_mismatch_count() == 1, "T2 counter");
    check(pose_error(motion.current_pose(), accepted) == 0.0, "T2 accepted pose unchanged");
    check(motion.current_time() == 3.0, "T2 accepted timestamp unchanged");
    check(!motion.raw_velocity_valid() && motion.body_twist().norm() == 0.0,
          "T2 velocity unchanged");
}

void t3_out_of_order()
{
    Motion motion;
    const auto twist = twist_cases()[1];
    const Eigen::Matrix4d origin = Eigen::Matrix4d::Identity();
    const Eigen::Matrix4d current = advanced(origin, twist, .010);
    motion.update(origin, 4.0, 2, 4.0);
    motion.update(current, 4.010, 2, 4.010);
    const auto velocity = motion.body_twist();
    check(motion.update(origin, 3.5, 2, 4.012) == Status::OutOfOrderTimestampRejected,
          "T3 status");
    check(motion.out_of_order_timestamp_count() == 1, "T3 counter");
    check(pose_error(motion.current_pose(), current) == 0.0, "T3 pose unchanged");
    check((motion.body_twist() - velocity).norm() == 0.0, "T3 velocity unchanged");
}

void t4_static_to_epoch()
{
    Motion motion;
    const auto twist = twist_cases()[2];
    const Eigen::Matrix4d staticPose = Eigen::Matrix4d::Identity();
    const Eigen::Matrix4d first = advanced(staticPose, twist, .010);
    const Eigen::Matrix4d second = advanced(first, twist, .010);
    motion.set_static_pose(staticPose);
    check(motion.update(first, 1788149018.000, 5, 1788149018.000)
              == Status::FirstTimestampedSample,
          "T4 first epoch sample seeds time");
    check(!motion.raw_velocity_valid() && motion.body_twist().norm() == 0.0,
          "T4 first epoch sample has no velocity");
    check(motion.last_elapsed() == 0.0, "T4 no epoch-zero elapsed");
    check(motion.update(second, 1788149018.010, 5, 1788149018.010)
              == Status::TimestampAdvanced,
          "T4 second source advances");
    check(std::abs(motion.last_elapsed() - .010) < 2e-8, "T4 real source interval");
    check((motion.body_twist() - twist).norm() < 2e-7, "T4 recovered epoch twist");
}

void t5_tiny_interval()
{
    Motion motion;
    const auto twist = twist_cases()[2];
    const Eigen::Matrix4d origin = Eigen::Matrix4d::Identity();
    const Eigen::Matrix4d tiny = advanced(origin, twist, .010);
    motion.update(origin, 5.0, 6, 5.0);
    check(motion.update(tiny, 5.0 + 1e-9, 6, 5.0 + 1e-9)
              == Status::TooSmallIntervalReseeded,
          "T5 tiny interval reseeded");
    check(motion.too_small_interval_count() == 1, "T5 counter");
    check(!motion.raw_velocity_valid() && motion.body_twist().norm() == 0.0,
          "T5 no huge twist");
    check(motion.current_pose().allFinite(), "T5 finite state");
    const Eigen::Matrix4d next = advanced(tiny, twist, .010);
    motion.update(next, 5.010000001, 6, 5.010000001);
    check((motion.body_twist() - twist).norm() < 2e-12,
          "T5 next valid interval reconstructs velocity");
}

void t6_freshness_and_t7_resume()
{
    Motion motion;
    const auto twist = twist_cases()[1];
    const Eigen::Matrix4d origin = Eigen::Matrix4d::Identity();
    const Eigen::Matrix4d measured = advanced(origin, twist, .010);
    motion.update(origin, 6.0, 7, 6.0);
    motion.update(measured, 6.010, 7, 6.010);
    for(const double age : {.049, .050})
    {
        motion.update(measured, 6.010, 7, 6.010 + age);
        check(motion.prediction_velocity_active(), "T6 velocity active at age <= 50 ms");
        check(pose_error(motion.predicted_pose(10, .002), advanced(measured, twist, .020)) < 2e-12,
              "T6 fresh pose predicts");
    }
    for(const double age : {.051, .100, .200})
    {
        motion.update(measured, 6.010, 7, 6.010 + age);
        check(!motion.prediction_velocity_active(), "T6 stale velocity inactive");
        check(motion.raw_velocity_valid(), "T6 raw velocity retained");
        check(pose_error(motion.predicted_pose(20, .002), measured) == 0.0,
              "T6 stale future pose held");
        check(predicted_parent_frame_state(motion, 20, .002).twistInBase.norm() == 0.0,
              "T6 stale predicted frame twist zero");
    }
    const Eigen::Matrix4d resumed = advanced(measured, twist, .210);
    check(motion.update(resumed, 6.220, 7, 6.220) == Status::TimestampAdvanced,
          "T7 source resumes");
    check((motion.body_twist() - twist).norm() < 2e-12,
          "T7 resumed velocity correct");
    check(motion.prediction_velocity_active(), "T7 prediction reactivated");
}

void t8_generation_reset()
{
    Motion motion;
    const auto twist = twist_cases()[0];
    const Eigen::Matrix4d a0 = Eigen::Matrix4d::Identity();
    const Eigen::Matrix4d a1 = advanced(a0, twist, .010);
    const Eigen::Matrix4d b0 = se3_exponential(
        (Eigen::Vector<double,6>() << .2, 0., 0., 0., 0., 0.).finished());
    const Eigen::Matrix4d b1 = advanced(b0, twist, .010);
    motion.update(a0, 100.000, 7, 100.000);
    motion.update(a1, 100.010, 7, 100.010);
    check(motion.update(b0, 0.000, 8, 0.000) == Status::GenerationResetFirstSample,
          "T8 generation-reset first status");
    check(!motion.raw_velocity_valid() && motion.body_twist().norm() == 0.0,
          "T8 no cross-generation velocity");
    check(motion.generation_reset_count() == 1 && motion.current_generation() == 8,
          "T8 generation observability");
    motion.update(b1, .010, 8, .010);
    check((motion.body_twist() - twist).norm() < 2e-12,
          "T8 generation-8 velocity");
}

void t9_nonfinite()
{
    Motion motion;
    const Eigen::Matrix4d origin = Eigen::Matrix4d::Identity();
    motion.update(origin, 1.0, 3, 1.0);
    for(int test = 0; test < 3; ++test)
    {
        Eigen::Matrix4d bad = origin;
        double timestamp = 1.010;
        double evaluation = 1.010;
        if(test == 0) bad(0,3) = std::numeric_limits<double>::quiet_NaN();
        if(test == 1) timestamp = std::numeric_limits<double>::infinity();
        if(test == 2) evaluation = std::numeric_limits<double>::quiet_NaN();
        bool threw = false;
        try { motion.update(bad, timestamp, 3, evaluation); }
        catch(const std::invalid_argument &) { threw = true; }
        check(threw, "T9 nonfinite throws");
        check(motion.last_update_status() == Status::NonfiniteInputRejected,
              "T9 nonfinite status");
        check(pose_error(motion.current_pose(), origin) == 0.0 && motion.current_time() == 1.0,
              "T9 accepted state unchanged");
    }
}

void t11_exact_stream_no_regression()
{
    Motion motion;
    const auto twist = twist_cases()[2];
    const Eigen::Matrix4d origin = Eigen::Matrix4d::Identity();
    for(int index = 0; index <= 50; ++index)
    {
        const double time = .002 * index;
        const Eigen::Matrix4d expected = advanced(origin, twist, time);
        motion.update(expected, time, 11, time);
        check(pose_error(motion.current_pose(), expected) < 2e-15,
              "T11 accepted pose exact");
        check(motion.measurement_age() == 0.0, "T11 zero age");
        if(index > 0)
        {
            check((motion.body_twist() - twist).norm() < 2e-12,
                  "T11 exact-stream twist");
            check(motion.prediction_velocity_active(), "T11 exact-stream active");
        }
    }
}

} // namespace

int main()
{
    t0_normal_constant_twist();
    t1_low_rate_source();
    t2_duplicate_pose_mismatch();
    t3_out_of_order();
    t4_static_to_epoch();
    t5_tiny_interval();
    t6_freshness_and_t7_resume();
    t8_generation_reset();
    t9_nonfinite();
    t11_exact_stream_no_regression();
    if(failures == 0)
    {
        std::cout << "parent_frame_temporal_semantics_test PASS\n";
        return 0;
    }
    std::cerr << "parent_frame_temporal_semantics_test FAIL (" << failures << ")\n";
    return 1;
}
