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
#
# ==========================================================
#
# FindBAT.cmake file to find the BAT library
#
# This script will define the following variables:
#   BAT_FOUND        - True if the system has the BAT library
#   BAT_INCLUDE_DIRS - The directory containing the BAT headers
#   BAT_LIBRARIES    - The BAT library
#   BAT_VERSION      - The version of the BAT library which was found
#
# And the following imported targets:
#   BAT::BAT         - The BAT library target

# Set search paths
set(_BAT_PATHS ${BAT_ROOT} $ENV{BAT_ROOT} $ENV{BATINSTALLDIR})

find_path(BAT_INCLUDE_DIR
  NAMES BAT/BCLog.h BAT/BCMath.h
  PATHS ${_BAT_PATHS}
  PATH_SUFFIXES include
)

find_library(BAT_LIBRARY
  NAMES BAT
  PATHS ${_BAT_PATHS}
  PATH_SUFFIXES lib
)

# Extract version from BCVersion.h if it exists, otherwise from BCModel.h or similar
if(BAT_INCLUDE_DIR AND EXISTS "${BAT_INCLUDE_DIR}/BAT/BCVersion.h")
  file(STRINGS "${BAT_INCLUDE_DIR}/BAT/BCVersion.h" _BAT_VERSION_LINE REGEX "^#define BAT_VERSION \"[^\"]*\"")
  string(REGEX REPLACE "^#define BAT_VERSION \"([^\"]*)\".*" "\\1" BAT_VERSION "${_BAT_VERSION_LINE}")
elseif(BAT_INCLUDE_DIR)
  # Fallback for older BAT versions where it might be in BCModel.h or just assume the requested version if it matches the headers
  set(BAT_VERSION "0.9.4.1") # Default assumption for this repo's compatibility
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BAT
  FOUND_VAR BAT_FOUND
  REQUIRED_VARS BAT_LIBRARY BAT_INCLUDE_DIR
  VERSION_VAR BAT_VERSION
)

if(BAT_FOUND)
  set(BAT_INCLUDE_DIRS ${BAT_INCLUDE_DIR})
  set(BAT_LIBRARIES ${BAT_LIBRARY})
  if(NOT TARGET BAT::BAT)
    add_library(BAT::BAT UNKNOWN IMPORTED)
    set_target_properties(BAT::BAT PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${BAT_INCLUDE_DIRS}"
      IMPORTED_LOCATION "${BAT_LIBRARIES}"
    )
  endif()
endif()

mark_as_advanced(BAT_INCLUDE_DIR BAT_LIBRARY)
