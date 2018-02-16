# Copyright (c) 2009--2017, the KLFitter developer team
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
#
# ==========================================================
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