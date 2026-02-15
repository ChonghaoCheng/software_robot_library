#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "RobotLibrary::Control" for configuration "Release"
set_property(TARGET RobotLibrary::Control APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(RobotLibrary::Control PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libControl.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS RobotLibrary::Control )
list(APPEND _IMPORT_CHECK_FILES_FOR_RobotLibrary::Control "${_IMPORT_PREFIX}/lib/libControl.a" )

# Import target "RobotLibrary::Math" for configuration "Release"
set_property(TARGET RobotLibrary::Math APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(RobotLibrary::Math PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libMath.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS RobotLibrary::Math )
list(APPEND _IMPORT_CHECK_FILES_FOR_RobotLibrary::Math "${_IMPORT_PREFIX}/lib/libMath.a" )

# Import target "RobotLibrary::Model" for configuration "Release"
set_property(TARGET RobotLibrary::Model APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(RobotLibrary::Model PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libModel.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS RobotLibrary::Model )
list(APPEND _IMPORT_CHECK_FILES_FOR_RobotLibrary::Model "${_IMPORT_PREFIX}/lib/libModel.a" )

# Import target "RobotLibrary::Trajectory" for configuration "Release"
set_property(TARGET RobotLibrary::Trajectory APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(RobotLibrary::Trajectory PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libTrajectory.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS RobotLibrary::Trajectory )
list(APPEND _IMPORT_CHECK_FILES_FOR_RobotLibrary::Trajectory "${_IMPORT_PREFIX}/lib/libTrajectory.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
