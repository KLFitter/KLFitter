# Set the minimum CMake version required to build the project.
cmake_minimum_required( VERSION 3.1 )

# Silence some warnings on macOS with new CMake versions.
if( NOT ${CMAKE_VERSION} VERSION_LESS 3.9 )
   cmake_policy( SET CMP0068 NEW )
endif()

# Set the project's name and version.
project( KLFitter )

# Enable using CTest in the project.
include( CTest )

# Set up the "C++ version" to use.
set( CMAKE_CXX_STANDARD_REQUIRED 11 CACHE STRING
   "Minimum C++ standard required for the build" )
set( CMAKE_CXX_STANDARD 14 CACHE STRING
   "C++ standard to use for the build" )
set( CMAKE_CXX_EXTENSIONS FALSE CACHE BOOL
   "(Dis)allow using compiler extensions" )

# If the user didn't request a build type explicitly, use an optimised
# build with debug symbols.
if( "${CMAKE_BUILD_TYPE}" STREQUAL "" )
   set( CMAKE_BUILD_TYPE "RelWithDebInfo" )
endif()
message( STATUS "Using build type: ${CMAKE_BUILD_TYPE}" )

# Set the warning flag(s) to use.
set( CMAKE_CXX_FLAGS "-Wall -pedantic" )

# Turn off the usage of RPATH completely:
set( CMAKE_SKIP_RPATH ON )
set( CMAKE_SKIP_BUILD_RPATH ON )
set( CMAKE_SKIP_INSTALL_RPATH ON )

# Set the location of the built libraries/executables inside the build
# directory.
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )

# Make the project's modules visible to CMake.
list( INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake )

# Add ROOT system directory and require ROOT.
find_package( ROOT REQUIRED )

# Figure out what to do with BAT.
set( _buildDir ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/BATBuild )
option( BUILTIN_BAT "Acquire BAT as part of building this project" OFF )
if( BUILTIN_BAT )
  # Download and install BAT using ExternalProject.
  include( ExternalProject )
  ExternalProject_Add( BAT
    PREFIX ${CMAKE_BINARY_DIR}/externals
    INSTALL_DIR ${CMAKE_BINARY_DIR}
    URL "https://github.com/bat/bat/releases/download/v0.9.4.1/BAT-0.9.4.1.tar.gz"
    URL_MD5 "e1a09dc5e379dd5b5d665729ea45b535"
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --with-rootsys=${ROOTSYS} --prefix=${_buildDir}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${_buildDir}/ <INSTALL_DIR> )
  # Set the include and library paths to use.
  set( BAT_ROOT $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}> )
  set( BAT_INCLUDE_DIR ${BAT_ROOT}/include )
  set( BAT_LIBRARY ${BAT_ROOT}/lib/libBAT.so )
  # Tell the user what will happen.
  message( STATUS "Using a privately downloaded BAT version for the build" )
else()
  # Look for BAT on the build system.
  find_package( BAT REQUIRED )
endif()

# Public header files for the shared/static library.
set( lib_headers include/KLFitter/BoostedLikelihoodTopLeptonJets.h
include/KLFitter/DetectorAtlas_7TeV.h include/KLFitter/DetectorAtlas_8TeV.h
include/KLFitter/DetectorBase.h include/KLFitter/Fitter.h
include/KLFitter/LikelihoodBase.h include/KLFitter/LikelihoodSgTopWtLJ.h
include/KLFitter/LikelihoodTTHLeptonJets.h
include/KLFitter/LikelihoodTTZTrilepton.h
include/KLFitter/LikelihoodTopAllHadronic.h
include/KLFitter/LikelihoodTopDilepton.h
include/KLFitter/LikelihoodTopLeptonJets.h
include/KLFitter/LikelihoodTopLeptonJetsUDSep.h
include/KLFitter/LikelihoodTopLeptonJets_Angular.h
include/KLFitter/LikelihoodTopLeptonJets_JetAngles.h
include/KLFitter/Particles.h include/KLFitter/Permutations.h )

# Source files for the shared/static library.
set( lib_sources src/BoostedLikelihoodTopLeptonJets.cxx
src/DetectorAtlas_7TeV.cxx src/DetectorAtlas_8TeV.cxx src/DetectorBase.cxx
src/Fitter.cxx src/LikelihoodBase.cxx src/LikelihoodSgTopWtLJ.cxx
src/LikelihoodTTHLeptonJets.cxx src/LikelihoodTTZTrilepton.cxx
src/LikelihoodTopAllHadronic.cxx src/LikelihoodTopDilepton.cxx
src/LikelihoodTopLeptonJets.cxx src/LikelihoodTopLeptonJetsUDSep.cxx
src/LikelihoodTopLeptonJets_Angular.cxx
src/LikelihoodTopLeptonJets_JetAngles.cxx src/Particles.cxx src/Permutations.cxx
src/PhysicsConstants.cxx src/ResDoubleGaussBase.cxx src/ResDoubleGaussE_1.cxx
src/ResDoubleGaussE_2.cxx src/ResDoubleGaussE_3.cxx src/ResDoubleGaussE_4.cxx
src/ResDoubleGaussE_5.cxx src/ResDoubleGaussPt.cxx src/ResGauss.cxx
src/ResGauss_MET.cxx src/ResolutionBase.cxx )

# Build the shared library.
add_library( KLFitter SHARED ${lib_headers} ${lib_sources} )
target_include_directories( KLFitter
   PUBLIC ${ROOT_INCLUDE_DIRS} ${BAT_INCLUDE_DIR}
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include> )
target_link_libraries( KLFitter ${ROOT_LIBRARIES} ${BAT_LIBARY} )
set_property( TARGET KLFitter
   PROPERTY PUBLIC_HEADER ${lib_headers} )
if( BUILTIN_BAT )
  add_dependencies( KLFitter BAT )
endif()

# Build the static library.
add_library( KLFitter-stat ${lib_headers} ${lib_sources} )
target_include_directories( KLFitter-stat
   PUBLIC ${ROOT_INCLUDE_DIRS} ${BAT_INCLUDE_DIR}
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include> )
target_link_libraries( KLFitter-stat ${ROOT_LIBRARIES} ${BAT_LIBRARY} )
set_property( TARGET KLFitter-stat
   PROPERTY PUBLIC_HEADER ${lib_headers} )
if( BUILTIN_BAT )
  add_dependencies( KLFitter-stat BAT )
endif()

# Install the libraries.
install( TARGETS KLFitter KLFitter-stat
   EXPORT KLFitterTargets
   ARCHIVE DESTINATION lib
   LIBRARY DESTINATION lib
   PUBLIC_HEADER DESTINATION include/KLFitter )

# Helper macro for building the project's unit tests.
macro( KLFitter_add_test name )
   # Build the unit-test executable:
   add_executable( ${name} ${ARGN} )
   target_link_libraries( ${name} KLFitter-stat )
   set_target_properties( ${name} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test-bin" )
   # Set up the test itself:
   add_test( NAME ${name}_ctest
      COMMAND ${CMAKE_BINARY_DIR}/test-bin/${name}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test-bin )
endmacro( KLFitter_add_test )

# Set up the test(s) of the project.
if( KLFITTER_INSTALL_TESTS )
  KLFitter_add_test( test-ljets-lh tests/test-ljets-lh.cxx )
endif()

# Install the CMake description of the project.
install( EXPORT KLFitterTargets
   FILE KLFitterConfig-targets.cmake
   DESTINATION cmake
   NAMESPACE "KLFitter::" )
