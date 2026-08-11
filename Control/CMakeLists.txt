
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
                    src/SerialLinkMovingFrameMPC.cpp
                    src/SerialLinkRMPCC.cpp
)

# Specify targets to be built
target_include_directories(Control PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/RobotLibrary>
)

target_link_libraries(Control PUBLIC Math Model Trajectory Eigen3::Eigen)                           # Other libraries needed to compile this one

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
