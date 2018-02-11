# FindBAT.cmake file to find the BAT library
#
# To ensure that this works correctly, please set the environment
# variable $BATINSTALLDIR to your BAT installation directory. This
# script will define the following variables:
#
#   BAT_FOUND    - True if the system has the BAT library
#   BAT_VERSION  - The version of the BAT library which was found
#
# and the following imported targets:
#
#   BAT::BAT   - The BAT library

find_path(BAT_INCLUDE_DIR
  NAMES BAT/BCLog.h BAT/BCMath.h
  PATHS ${CMAKE_INSTALL_PREFIX}/include $ENV{BATINSTALLDIR}/include
)
find_library(BAT_LIBRARY
  NAMES BAT
  PATHS ${CMAKE_INSTALL_PREFIX}/lib $ENV{BATINSTALLDIR}/lib
)

set(BAT_VERSION ${PC_BAT_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BAT
  FOUND_VAR BAT_FOUND
  REQUIRED_VARS
    BAT_LIBRARY
    BAT_INCLUDE_DIR
  VERSION_VAR BAT_VERSION
)