# Copyright (c) 2009--2018, the KLFitter developer team
#
# This file is part of KLFitter.
#
# KLFitter is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# KLFitter is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with KLFitter. If not, see <http://www.gnu.org/licenses/>.

message( STATUS "Installing KLFitter test suite ..." )

# Check whether a privately downloaded version of googletest was requested.
if( KLF_BUILD_GTEST )
  include( ExternalProject )
  set( _buildDir ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/GTestBuild )
  ExternalProject_Add( googletest_ext
    PREFIX ${CMAKE_BINARY_DIR}/externals
    INSTALL_DIR ${CMAKE_BINARY_DIR}
    URL "https://github.com/google/googletest/archive/release-1.8.1.tar.gz"
    URL_HASH SHA256=9bf1fe5182a604b4135edc1a425ae356c9ad15e9b23f9f12a02e80184c3a249c
    CMAKE_ARGS
    -DBUILD_GMOCK:BOOL=Off
    -DGTEST_HAS_PTHREAD=0
    -DCMAKE_INSTALL_PREFIX:PATH=${_buildDir}
  )

  # Add an additional step to copy googletest into the install directory.
  ExternalProject_Add_Step( googletest_ext CopyInstall
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${_buildDir}/ <INSTALL_DIR>
    DEPENDEES install
  )

  # Add the dependency on the thread library (pthread).
  find_package( Threads )

  # Set the include directory for googletest.
  ExternalProject_Get_Property( googletest_ext install_dir )
  set( GTEST_INCLUDE_DIR ${install_dir}/include )
  set( GTEST_INCLUDE_DIRS ${install_dir}/include )

  # Set the library path and library name for googletest.
  get_property( _use_lib64_paths GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS )
  if( _use_lib64_paths )
    set( _lib_path "${install_dir}/lib64" )
  else()
    set( _lib_path "${install_dir}/lib" )
  endif()
  set( GTEST_LIBRARIES "${_lib_path}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a" )
  set( GTEST_MAIN_LIBRARIES "${_lib_path}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main.a" )
  set( GTEST_BOTH_LIBRARIES "${GTEST_LIBRARIES};${GTEST_MAIN_LIBRARIES}" )

  # Now add the GTest and GTestMain libraries including their location.
  add_library( GTest::GTest UNKNOWN IMPORTED )
  set_property( TARGET GTest::GTest PROPERTY IMPORTED_LOCATION "${GTEST_LIBRARIES}" )
  add_library( GTest::GTestMain UNKNOWN IMPORTED)
  set_property( TARGET GTest::GTestMain PROPERTY IMPORTED_LOCATION "${GTEST_MAIN_LIBRARIES}" )

  add_dependencies( GTest::GTest googletest_ext )
  add_dependencies( GTest::GTest ${CMAKE_THREAD_LIBS_INIT} )

  # Write status message for the user.
  message( STATUS "Using a privately downloaded googletest version for the build" )
else()
  find_package( GTest 1.8.1 REQUIRED )
endif()

# Build the main unit-test executable.
add_executable( gtests.exe gtests.cxx )
target_include_directories( gtests.exe
  PUBLIC ${ROOT_INCLUDE_DIRS} ${BAT_INCLUDE_DIR} ${GTEST_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/tests>
)
target_link_libraries( gtests.exe KLFitter GTest::GTest ${CMAKE_THREAD_LIBS_INIT})
set_target_properties( gtests.exe PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test-bin" )
if( KLF_BUILD_GTEST )
  add_dependencies( gtests.exe googletest_ext )
endif()

# Build the integration test for the l+jets likelihood.
add_executable( test-ljets-lh.exe test-ljets-lh.cxx )
target_link_libraries( test-ljets-lh.exe KLFitter )
set_target_properties( test-ljets-lh.exe PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test-bin" )
