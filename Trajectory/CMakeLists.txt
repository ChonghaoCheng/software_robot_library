# Trajectory/CMakeLists.txt
add_library(Trajectory src/CartesianSpline.cpp
                       src/MinimumArcLength.cpp
                       src/SplineTrajectory.cpp
                       src/TrajectoryBase.cpp
                       src/TrapezoidalVelocity.cpp
)

target_include_directories(Trajectory PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/RobotLibrary>
)

target_link_libraries(Trajectory PUBLIC Math Model Eigen3::Eigen)

if(BUILD_TESTING)
    add_executable(cartesian_trajectory_frame_test test/CartesianTrajectoryFrame.cpp)
    target_link_libraries(cartesian_trajectory_frame_test PRIVATE Trajectory)
    add_test(NAME cartesian_trajectory_frame_test COMMAND cartesian_trajectory_frame_test)

    add_executable(cartesian_spline_quaternion_hemisphere_test
                   test/CartesianSplineQuaternionHemisphere.cpp)
    target_link_libraries(cartesian_spline_quaternion_hemisphere_test PRIVATE Trajectory)
    add_test(NAME cartesian_spline_quaternion_hemisphere_test
             COMMAND cartesian_spline_quaternion_hemisphere_test)

    add_executable(cartesian_spline_pose_only_test
                   test/CartesianSplinePoseOnly.cpp)
    target_link_libraries(cartesian_spline_pose_only_test PRIVATE Trajectory)
    add_test(NAME cartesian_spline_pose_only_test
             COMMAND cartesian_spline_pose_only_test)
endif()

install(TARGETS  Trajectory
        EXPORT   TrajectoryTargets
        LIBRARY  DESTINATION lib
        ARCHIVE  DESTINATION lib
        RUNTIME  DESTINATION bin
        INCLUDES DESTINATION include
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
        DESTINATION include/RobotLibrary
        FILES_MATCHING PATTERN "*.h"
)
