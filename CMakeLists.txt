# Auto-generated on: 2017-03-09 17:44:32.719752

# Declare the name of this package:
atlas_subdir( KLFitter None )

# This package uses ROOT:
#find_package( ROOT REQUIRED COMPONENTS Core Gpad Tree Hist RIO Math MathMore MathCore Graf Physics )
find_package( ROOT REQUIRED  COMPONENTS Core RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread Minuit )

MESSAGE( STATUS "CMAKE_MODULE_PATH : " ${CMAKE_MODULE_PATH} )
find_package( BAT REQUIRED )
MESSAGE( STATUS "BAT_LIBRARIES used in KLFitter : " ${BAT_LIBRARIES} )

# Build a library that other components can link against:
atlas_add_library( KLFitter Root/*.cxx Root/*.h Root/*.icc
                   KLFitter/*.h KLFitter/*.icc KLFitter/*/*.h
                   KLFitter/*/*.icc 
                   PUBLIC_HEADERS KLFitter
                   LINK_LIBRARIES ${ROOT_LIBRARIES}
				  ${BAT_LIBRARIES}
                   INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
		     		${BAT_INCLUDE_DIRS})

# Install data files from the package:
atlas_install_data( data/* )
