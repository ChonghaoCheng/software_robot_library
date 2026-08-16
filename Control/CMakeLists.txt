
# List the source files for this library
add_library(Control src/AdmittanceContactController.cpp
                    src/ContactControllerBase.cpp
                    src/SerialLinkBase.cpp
                    src/SerialLinkGeometricImpedance.cpp
                    src/SerialLinkImpedance.cpp
                    src/SerialLinkKinematic.cpp
                    src/SerialLinkVelocityBase.cpp
                    src/SerialLinkMPCC.cpp
                    src/SerialLinkCartesianMPCC.cpp
                    src/SerialLinkMPC.cpp
                    src/SerialLinkLieAlgebraMPC.cpp
                    src/SerialLinkMovingFrameMPC.cpp
                    src/SerialLinkRMPCC.cpp
)

# Specify targets to be built
target_include_directories(Control PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/RobotLibrary>
)

target_link_libraries(Control PUBLIC Math Model Trajectory Eigen3::Eigen)                           # Other libraries needed to compile this one

if(BUILD_TESTING)
    add_executable(rmpcc_prediction_test test/RmpccPrediction.cpp)
    target_link_libraries(rmpcc_prediction_test PRIVATE Control)
    add_test(NAME rmpcc_prediction_test COMMAND rmpcc_prediction_test)

    add_executable(lie_algebra_mpc_prediction_test test/LieAlgebraMPCPrediction.cpp)
    target_link_libraries(lie_algebra_mpc_prediction_test PRIVATE Control)
    add_test(NAME lie_algebra_mpc_prediction_test COMMAND lie_algebra_mpc_prediction_test)

    add_executable(cartesian_mpcc_prediction_test test/CartesianMPCCPrediction.cpp)
    target_link_libraries(cartesian_mpcc_prediction_test PRIVATE Control)
    add_test(NAME cartesian_mpcc_prediction_test COMMAND cartesian_mpcc_prediction_test)

    add_executable(rmpcc_cost_geometry_test test/RmpccCostGeometry.cpp)
    target_link_libraries(rmpcc_cost_geometry_test PRIVATE Control)
    add_test(NAME rmpcc_cost_geometry_test COMMAND rmpcc_cost_geometry_test)

    add_executable(rmpcc_full_residual_linearization_test
                   test/RmpccFullResidualLinearization.cpp)
    target_link_libraries(rmpcc_full_residual_linearization_test PRIVATE Control)
    add_test(NAME rmpcc_full_residual_linearization_test
             COMMAND rmpcc_full_residual_linearization_test)

    add_executable(rmpcc_full_objective_directional_derivative_test
                   test/RmpccFullObjectiveDirectionalDerivative.cpp)
    target_link_libraries(rmpcc_full_objective_directional_derivative_test
                          PRIVATE Control)
    add_test(NAME rmpcc_full_objective_directional_derivative_test
             COMMAND rmpcc_full_objective_directional_derivative_test)

    add_executable(rmpcc_decoupled_residual_linearization_test
                   test/RmpccDecoupledResidualLinearization.cpp)
    target_link_libraries(rmpcc_decoupled_residual_linearization_test PRIVATE Control)
    add_test(NAME rmpcc_decoupled_residual_linearization_test
             COMMAND rmpcc_decoupled_residual_linearization_test)

    add_executable(rmpcc_schedule_limit_test test/RmpccScheduleLimit.cpp)
    target_link_libraries(rmpcc_schedule_limit_test PRIVATE Control)
    add_test(NAME rmpcc_schedule_limit_test COMMAND rmpcc_schedule_limit_test)

    add_executable(mpcc_stage_tangent_coordinates_test
                   test/MpccStageTangentCoordinates.cpp)
    target_link_libraries(mpcc_stage_tangent_coordinates_test PRIVATE Control)
    add_test(NAME mpcc_stage_tangent_coordinates_test
             COMMAND mpcc_stage_tangent_coordinates_test)

    add_executable(common_progress_schedule_test test/CommonProgressSchedule.cpp)
    target_link_libraries(common_progress_schedule_test PRIVATE Control)
    add_test(NAME common_progress_schedule_test COMMAND common_progress_schedule_test)
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
