
# Keep one binary target while organizing its implementation by controller purpose.
add_library(Control src/Core/SerialLinkBase.cpp
                    src/Core/SerialLinkVelocityBase.cpp
                    src/TrajectoryTracking/SerialLinkKinematic.cpp
                    src/TrajectoryTracking/SerialLinkMPCC.cpp
                    src/TrajectoryTracking/SerialLinkCartesianMPCC.cpp
                    src/TrajectoryTracking/SerialLinkTimeIndexedMPC.cpp
                    src/TrajectoryTracking/SerialLinkMPC.cpp
                    src/TrajectoryTracking/SerialLinkLieAlgebraMPC.cpp
                    src/TrajectoryTracking/SerialLinkRMPCC.cpp
                    src/Contact/AdmittanceContactController.cpp
                    src/Contact/ContactControllerBase.cpp
                    src/Contact/SerialLinkGeometricImpedance.cpp
                    src/Contact/SerialLinkImpedance.cpp
                    src/Contact/SerialLinkMovingFrameMPC.cpp
)

# Specify targets to be built
target_include_directories(Control PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/RobotLibrary>
)

target_link_libraries(Control PUBLIC Math Model Trajectory Eigen3::Eigen)                           # Other libraries needed to compile this one

if(BUILD_TESTING)
    function(add_control_test target source)
        add_executable(${target} ${source})
        target_link_libraries(${target} PRIVATE Control)
        target_include_directories(${target} PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src/TrajectoryTracking/detail)
        add_test(NAME ${target} COMMAND ${target})
    endfunction()

    add_control_test(rmpcc_prediction_test test/TrajectoryTracking/RmpccPrediction.cpp)
    add_control_test(rmpcc_reference_motion_test test/TrajectoryTracking/RmpccReferenceMotion.cpp)
    add_control_test(lie_algebra_mpc_prediction_test test/TrajectoryTracking/LieAlgebraMPCPrediction.cpp)
    add_control_test(time_indexed_mpc_diagnostics_test test/TrajectoryTracking/TimeIndexedMpcDiagnostics.cpp)
    add_control_test(cartesian_mpcc_prediction_test test/TrajectoryTracking/CartesianMPCCPrediction.cpp)
    add_control_test(rmpcc_cost_geometry_test test/TrajectoryTracking/RmpccCostGeometry.cpp)
    add_control_test(rmpcc_full_residual_linearization_test test/TrajectoryTracking/RmpccFullResidualLinearization.cpp)
    add_control_test(rmpcc_full_objective_directional_derivative_test test/TrajectoryTracking/RmpccFullObjectiveDirectionalDerivative.cpp)
    add_control_test(rmpcc_phase_residual_test test/TrajectoryTracking/RmpccPhaseResidual.cpp)
    add_control_test(rmpcc_residual_non_expansion_test test/TrajectoryTracking/RmpccResidualNonExpansion.cpp)
    add_control_test(rmpcc_associated_phase_conditioning_test test/TrajectoryTracking/RmpccAssociatedPhaseConditioning.cpp)
    add_control_test(rmpcc_phase_association_profiles_test test/TrajectoryTracking/RmpccPhaseAssociationProfiles.cpp)
    add_control_test(n126_associated_phase_test test/TrajectoryTracking/N126AssociatedPhase.cpp)
    target_include_directories(n126_associated_phase_test PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/TrajectoryTracking/detail)
    add_control_test(rmpcc_decoupled_residual_linearization_test test/TrajectoryTracking/RmpccDecoupledResidualLinearization.cpp)
    add_control_test(rmpcc_schedule_limit_test test/TrajectoryTracking/RmpccScheduleLimit.cpp)
    add_control_test(mpcc_stage_tangent_coordinates_test test/TrajectoryTracking/MpccStageTangentCoordinates.cpp)
    add_control_test(common_progress_schedule_test test/TrajectoryTracking/CommonProgressSchedule.cpp)
    add_control_test(rmpcc_qp_constraints_test test/TrajectoryTracking/RmpccQpConstraints.cpp)

    # Offline N125-QP diagnostic utility. It is intentionally not a CTest
    # because it requires an explicitly supplied, non-committed snapshot.
    add_executable(n125_qp_active_set_replay
        test/TrajectoryTracking/N125QpActiveSetReplay.cpp)
    target_link_libraries(n125_qp_active_set_replay PRIVATE Control)

    # Deterministic theory-to-mechanism study. Formal execution supplies an
    # output directory; keeping it out of CTest avoids regenerating evidence
    # during ordinary unit-test runs.
    add_executable(n125_task_phase_validation
        test/TrajectoryTracking/N125TaskPhaseValidation.cpp)
    target_link_libraries(n125_task_phase_validation PRIVATE Control)
    target_include_directories(n125_task_phase_validation PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/TrajectoryTracking/detail)

    # DIAG-01 exact-on-path identity gate. Formal execution supplies the two
    # frozen task-frame path CSVs and an evidence output path.
    add_executable(n135_diag01_zero_point
        test/TrajectoryTracking/N135Diag01ZeroPoint.cpp)
    target_link_libraries(n135_diag01_zero_point PRIVATE Control)
    target_include_directories(n135_diag01_zero_point PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/TrajectoryTracking/detail)

    # Offline solver-component timing utility used by PERF-01. It is not a
    # CTest because timings are machine/environment dependent.
    add_executable(rmpcc_perf01_benchmark
        test/TrajectoryTracking/RmpccPerf01Benchmark.cpp)
    target_link_libraries(rmpcc_perf01_benchmark PRIVATE Control)
endif()

# Installation instructions
install(TARGETS  Control
        EXPORT   ControlTargets
        LIBRARY  DESTINATION lib
        ARCHIVE  DESTINATION lib
        RUNTIME  DESTINATION bin
        INCLUDES DESTINATION include
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/                                              # The directory we want to install
        DESTINATION include/RobotLibrary                                                            # Where to install it (no trailing slash, matches other components)
        FILES_MATCHING PATTERN "*.h"
)
