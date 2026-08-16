
# Keep one binary target while organizing its implementation by controller purpose.
add_library(Control src/Core/SerialLinkBase.cpp
                    src/Core/SerialLinkVelocityBase.cpp
                    src/TrajectoryTracking/SerialLinkKinematic.cpp
                    src/TrajectoryTracking/SerialLinkMPCC.cpp
                    src/TrajectoryTracking/SerialLinkCartesianMPCC.cpp
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
    add_control_test(lie_algebra_mpc_prediction_test test/TrajectoryTracking/LieAlgebraMPCPrediction.cpp)
    add_control_test(cartesian_mpcc_prediction_test test/TrajectoryTracking/CartesianMPCCPrediction.cpp)
    add_control_test(rmpcc_cost_geometry_test test/TrajectoryTracking/RmpccCostGeometry.cpp)
    add_control_test(rmpcc_full_residual_linearization_test test/TrajectoryTracking/RmpccFullResidualLinearization.cpp)
    add_control_test(rmpcc_full_objective_directional_derivative_test test/TrajectoryTracking/RmpccFullObjectiveDirectionalDerivative.cpp)
    add_control_test(rmpcc_phase_residual_test test/TrajectoryTracking/RmpccPhaseResidual.cpp)
    add_control_test(rmpcc_decoupled_residual_linearization_test test/TrajectoryTracking/RmpccDecoupledResidualLinearization.cpp)
    add_control_test(rmpcc_schedule_limit_test test/TrajectoryTracking/RmpccScheduleLimit.cpp)
    add_control_test(mpcc_stage_tangent_coordinates_test test/TrajectoryTracking/MpccStageTangentCoordinates.cpp)
    add_control_test(common_progress_schedule_test test/TrajectoryTracking/CommonProgressSchedule.cpp)
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
