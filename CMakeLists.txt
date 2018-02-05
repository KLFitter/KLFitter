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
set( CMAKE_CXX_STANDARD 11 CACHE STRING
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
# list( INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake )

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
   PUBLIC ${EIGEN3_INCLUDE_DIR}
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>
   PRIVATE ${Boost_INCLUDE_DIRS} )
set_property( TARGET KLFitter
   PROPERTY PUBLIC_HEADER ${lib_headers} )

# Build the static library.
add_library( KLFitter-stat ${lib_headers} ${lib_sources} )
target_include_directories( KLFitter-stat
   PUBLIC ${EIGEN3_INCLUDE_DIR}
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>
   PRIVATE ${Boost_INCLUDE_DIRS} )
set_property( TARGET KLFitter-stat
   PROPERTY PUBLIC_HEADER ${lib_headers} )

# Install the libraries.
install( TARGETS KLFitter KLFitter-stat
   EXPORT KLFitterTargets
   ARCHIVE DESTINATION lib
   LIBRARY DESTINATION lib
   PUBLIC_HEADER DESTINATION include/KLFitter )

# Install the CMake description of the project.
install( EXPORT KLFitterTargets
   FILE KLFitterConfig-targets.cmake
   DESTINATION cmake
   NAMESPACE "KLFitter::" )
