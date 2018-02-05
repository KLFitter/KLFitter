#
# Configuration file allowing clients to build code against the KLFitter
# library/libraries.
#

# Figure out the root directory of the installation.
get_filename_component( _thisdir "${CMAKE_CURRENT_LIST_FILE}" PATH )
get_filename_component( _basedir "${_thisdir}" PATH )
get_filename_component( KLFitter_INSTALL_DIR "${_basedir}" ABSOLUTE CACHE )

# Tell the user what happened.
if( NOT KLFitter_FIND_QUIETLY )
   message( STATUS
      "Found KLFitter: ${KLFitter_INSTALL_DIR} (version: ${KLFitter_VERSION})" )
endif()

# Include the "targets file".
# include( ${_thisdir}/KLFitterConfig-targets.cmake )

# Set some "old style" variables for using the KLFitter installation.
set( KLFitter_INCLUDE_DIRS "${KLFitter_INSTALL_DIR}/include" )
set( KLFitter_LIBRARIES KLFitter::KLFitter )

# Clean up.
unset( _thisdir )
unset( _basedir )