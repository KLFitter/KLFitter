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

# Check whether a privately downloaded version of googletest was requested.
if( KLF_BUILD_GTEST )
  include( ExternalProject )
  set( _buildDir ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/GTestBuild )
  ExternalProject_Add( googletest_ext
    PREFIX ${CMAKE_BINARY_DIR}/externals
    INSTALL_DIR ${CMAKE_BINARY_DIR}
    URL "https://github.com/google/googletest/archive/refs/tags/release-1.10.0.tar.gz"
    URL_HASH SHA256=9dc9157a9a1551ec7a7e43daea9a694a0bb5fb8bec81235d8a1e6ef64c716dcb
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
  find_package( GTest 1.10.0 REQUIRED )
endif()

# List of all unit-test implementation files.
set( test_files
  units/Particles/TestBoson.cxx
  units/Particles/TestElectron.cxx
  units/Particles/TestMuon.cxx
  units/Particles/TestNeutrino.cxx
  units/Particles/TestParton.cxx
  units/Particles/TestPhoton.cxx
  units/Particles/TestTau.cxx
  units/Particles/TestTrack.cxx
  units/TestParticleCollection.cxx
)

# Build the main unit-test executable.
add_executable( unit_tests.exe unit_tests.cxx ${test_files} )
target_include_directories( unit_tests.exe
  PUBLIC ${ROOT_INCLUDE_DIRS} ${GTEST_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_link_libraries( unit_tests.exe GTest::GTest KLFitter ${CMAKE_THREAD_LIBS_INIT})
set_target_properties( unit_tests.exe PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test-bin" )
if( KLF_BUILD_GTEST )
  add_dependencies( unit_tests.exe googletest_ext )
endif()

# Build the integration test for the l+jets likelihood.
add_executable( test-ljets-lh.exe integration/test-ljets-lh.cxx )
target_link_libraries( test-ljets-lh.exe KLFitter )
set_target_properties( test-ljets-lh.exe PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test-bin" )
