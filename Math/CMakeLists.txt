# Math/CMakeLists.txt
add_library(Math src/CondensedMPC.cpp
                 src/MathFunctions.cpp
                 src/Polynomial.cpp
                 src/SkewSymmetric.cpp
                 src/Spline.cpp
)

target_include_directories(Math PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/RobotLibrary>
)

target_link_libraries(Math PUBLIC Eigen3::Eigen)

if(BUILD_TESTING)
    add_executable(qp_active_set_test test/QPSolverActiveSet.cpp)
    target_link_libraries(qp_active_set_test PRIVATE Math)
    add_test(NAME qp_active_set_test COMMAND qp_active_set_test)

    add_executable(se3_right_jacobian_test test/SE3RightJacobian.cpp)
    target_link_libraries(se3_right_jacobian_test PRIVATE Math)
    add_test(NAME se3_right_jacobian_test COMMAND se3_right_jacobian_test)

    add_executable(condensed_time_varying_mpc_test test/CondensedTimeVaryingMPC.cpp)
    target_link_libraries(condensed_time_varying_mpc_test PRIVATE Math)
    add_test(NAME condensed_time_varying_mpc_test COMMAND condensed_time_varying_mpc_test)

    add_executable(discrete_integrator_lqr_test test/DiscreteIntegratorLQR.cpp)
    target_link_libraries(discrete_integrator_lqr_test PRIVATE Math)
    add_test(NAME discrete_integrator_lqr_test COMMAND discrete_integrator_lqr_test)
endif()

install(TARGETS  Math
        EXPORT   MathTargets
        LIBRARY  DESTINATION lib
        ARCHIVE  DESTINATION lib
        RUNTIME  DESTINATION bin
        INCLUDES DESTINATION include
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
        DESTINATION include/RobotLibrary
        FILES_MATCHING PATTERN "*.h"
)
