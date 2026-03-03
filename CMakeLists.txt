# Copyright (c) 2009--2026, the KLFitter developer team
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

# Set the minimum CMake version required to build the project.
cmake_minimum_required( VERSION 3.10 )

# Set the project's name and version.
project( KLFitter VERSION 1.3.2 )

# Enable using CTest in the project.
include( CTest )

# Set up the "C++ version" to use.
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_STANDARD 14 )
set( CMAKE_CXX_EXTENSIONS FALSE )

# Default to building shared libraries.
if( NOT DEFINED BUILD_SHARED_LIBS )
   set( BUILD_SHARED_LIBS ON )
endif()

# If the user didn't request a build type explicitly, use an optimised
# build with debug symbols.
if( "${CMAKE_BUILD_TYPE}" STREQUAL "" )
   set( CMAKE_BUILD_TYPE "RelWithDebInfo" )
endif()
message( STATUS "Using build type: ${CMAKE_BUILD_TYPE}" )

# Set the warning flag(s) to use.
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic" )

# Set up RPATH so that executables find the library automatically.
set( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
set( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )
list( FIND CMAKE_PLATFORM_STAGED_RUNTIME_LIBS "${CMAKE_INSTALL_PREFIX}/lib" _isSystem )
if( _isSystem STREQUAL "-1" )
   set( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
endif()

# Set the location of the built libraries/executables inside the build
# directory.
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )

# Make the project's modules visible to CMake.
list( INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake )

# Add ROOT system directory and require ROOT.
find_package( ROOT 5.34.10 REQUIRED MathCore Minuit )

# Figure out what to do with BAT.
set( _buildDir ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/BATBuild )
option( BUILTIN_BAT "Acquire BAT as part of building this project" OFF )
message( STATUS "BUILTIN_BAT is set to: ${BUILTIN_BAT}" )
if( BUILTIN_BAT )
  # Download and install BAT using ExternalProject.
  include( ExternalProject )
  ExternalProject_Add( BAT
    PREFIX ${CMAKE_BINARY_DIR}/externals
    URL "https://github.com/bat/bat/releases/download/v0.9.4.1/BAT-0.9.4.1.tar.gz"
    URL_HASH SHA256=D46C6F834CB5888BBF4DB393887190380132FA48816E0804F79C4A3CC344EF87
    PATCH_COMMAND /bin/bash ${CMAKE_SOURCE_DIR}/external/bat/patch_bat.sh <SOURCE_DIR> ${CMAKE_SOURCE_DIR}/external/bat
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --with-rootsys=${ROOTSYS} --prefix=${CMAKE_INSTALL_PREFIX} --enable-silent-rules
    BUILD_COMMAND make LIBTOOLFLAGS=--silent
    INSTALL_COMMAND make install LIBTOOLFLAGS=--silent )
  # Set the include and library paths to use for the build.
  set( BAT_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include )
  set( BAT_LIBRARY ${CMAKE_INSTALL_PREFIX}/lib/libBAT.so )
  # Tell the user what will happen.
  message( STATUS "Using a privately downloaded BAT version for the build" )
else()
  # Look for BAT on the build system.
  find_package( BAT 0.9.4.1 EXACT REQUIRED )
  message( STATUS "Found BAT: ${BAT_LIBRARIES} (version ${BAT_VERSION})" )
endif()

# Public header files for the shared/static library.
set( lib_headers
  include/KLFitter/BoostedLikelihoodTopLeptonJets.h
  include/KLFitter/DetectorAtlas_7TeV.h
  include/KLFitter/DetectorAtlas_8TeV.h
  include/KLFitter/DetectorAtlas_CrystalBall.h
  include/KLFitter/DetectorSnowmass.h
  include/KLFitter/DetectorBase.h
  include/KLFitter/Fitter.h
  include/KLFitter/LikelihoodBase.h
  include/KLFitter/LikelihoodOneHadronicTop.h
  include/KLFitter/LikelihoodSgTopWtLJ.h
  include/KLFitter/LikelihoodTTHLeptonJets.h
  include/KLFitter/LikelihoodTTZTrilepton.h
  include/KLFitter/LikelihoodTopAllHadronic.h
  include/KLFitter/LikelihoodTopDilepton.h
  include/KLFitter/LikelihoodTopLeptonJets.h
  include/KLFitter/LikelihoodTopLeptonJetsUDSep.h
  include/KLFitter/LikelihoodTopLeptonJets_Angular.h
  include/KLFitter/LikelihoodTopLeptonJets_JetAngles.h
  include/KLFitter/LikelihoodTwoTracks.h
  include/KLFitter/Particles.h
  include/KLFitter/Permutations.h
  include/KLFitter/PhysicsConstants.h )

# Source files for the shared/static library.
set( lib_sources
  src/BoostedLikelihoodTopLeptonJets.cxx
  src/DetectorAtlas_7TeV.cxx
  src/DetectorAtlas_8TeV.cxx
  src/DetectorAtlas_CrystalBall.cxx
  src/DetectorSnowmass.cxx
  src/DetectorBase.cxx
  src/Fitter.cxx
  src/LikelihoodBase.cxx
  src/LikelihoodOneHadronicTop.cxx
  src/LikelihoodSgTopWtLJ.cxx
  src/LikelihoodTTHLeptonJets.cxx
  src/LikelihoodTTZTrilepton.cxx
  src/LikelihoodTopAllHadronic.cxx
  src/LikelihoodTopDilepton.cxx
  src/LikelihoodTopLeptonJets.cxx
  src/LikelihoodTopLeptonJetsUDSep.cxx
  src/LikelihoodTopLeptonJets_Angular.cxx
  src/LikelihoodTopLeptonJets_JetAngles.cxx
  src/LikelihoodTwoTracks.cxx
  src/Particles.cxx
  src/Permutations.cxx
  src/PhysicsConstants.cxx
  src/ResCrystalBallBase.cxx
  src/ResCrystalBallJets.cxx
  src/ResDoubleGaussBase.cxx
  src/ResDoubleGaussE_1.cxx
  src/ResDoubleGaussE_2.cxx
  src/ResDoubleGaussE_3.cxx
  src/ResDoubleGaussE_4.cxx
  src/ResDoubleGaussE_5.cxx
  src/ResDoubleGaussPt.cxx
  src/ResGauss.cxx
  src/ResGaussE.cxx
  src/ResGaussPt.cxx
  src/ResGauss_MET.cxx
  src/ResSingleGaussBase.cxx
  src/ResSingleGaussLinearBase.cxx
  src/ResSingleGaussE.cxx
  src/ResSingleGaussMET.cxx
  src/ResSingleGaussPt.cxx
  src/ResolutionBase.cxx )

# Build the library.
add_library( KLFitter ${lib_headers} ${lib_sources} )
add_library( KLFitter::KLFitter ALIAS KLFitter )
target_compile_features( KLFitter PUBLIC cxx_std_14 )
target_include_directories( KLFitter
   PUBLIC ${ROOT_INCLUDE_DIRS} ${BAT_INCLUDE_DIR}
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include> )
target_link_libraries( KLFitter ${ROOT_LIBRARIES} ${BAT_LIBRARY} )
set_property( TARGET KLFitter
   PROPERTY PUBLIC_HEADER ${lib_headers} )
if( BUILTIN_BAT )
  add_dependencies( KLFitter BAT )
endif()

# Install the library.
install( TARGETS KLFitter
   EXPORT KLFitterTargets
   ARCHIVE DESTINATION lib
   LIBRARY DESTINATION lib
   PUBLIC_HEADER DESTINATION include/KLFitter )

# Install the transfer functions.
install( DIRECTORY data/transferfunctions/
  DESTINATION share/KLFitter/transferfunctions )

# Helper macro for building the project's unit tests.
macro( KLFitter_add_test name )
   # Build the unit-test executable:
   add_executable( ${name} ${ARGN} )
   target_link_libraries( ${name} KLFitter )
   set_target_properties( ${name} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test-bin" )
   # Set up the test itself:
   add_test( NAME ${name}_ctest
      COMMAND ${CMAKE_BINARY_DIR}/test-bin/${name}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test-bin )
endmacro( KLFitter_add_test )

# Set up the test(s) of the project.
option( INSTALL_TESTS "Install the unit tests to validate KLFitter installation" OFF )
if( INSTALL_TESTS )
  KLFitter_add_test( test-ljets-lh.exe tests/test-ljets-lh.cxx )
endif()

# Helper macro for building the project's executables.
macro( KLFitter_add_executable name )
  add_executable( ${name} ${ARGN} )
  target_link_libraries( ${name} KLFitter )
  install( TARGETS ${name}
    EXPORT KLFitterTargets
    RUNTIME DESTINATION bin )
endmacro( KLFitter_add_executable )

# Build the example executable(s) of the project.
option( INSTALL_EXAMPLES "Install the example executable(s) of KLFitter" ON )
if( INSTALL_EXAMPLES )
  KLFitter_add_executable( example-top-ljets.exe util/example-top-ljets.cxx )
endif()

# Install the CMake description of the project.
install( EXPORT KLFitterTargets
   FILE KLFitterTargets.cmake
   DESTINATION lib/cmake/KLFitter
   NAMESPACE "KLFitter::" )

# Create a version file for the project.
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/KLFitterConfigVersion.cmake.in"
   "${CMAKE_CURRENT_BINARY_DIR}/KLFitterConfigVersion.cmake" @ONLY )

# Create and install the main KLFitterConfig.cmake file.
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/KLFitterConfig.cmake.in"
   "${CMAKE_CURRENT_BINARY_DIR}/KLFitterConfig.cmake" @ONLY )

install( FILES
   "${CMAKE_CURRENT_BINARY_DIR}/KLFitterConfig.cmake"
   "${CMAKE_CURRENT_BINARY_DIR}/KLFitterConfigVersion.cmake"
   DESTINATION lib/cmake/KLFitter )
