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
  PATHS ${BAT_ROOT} $ENV{BATINSTALLDIR}
  PATH_SUFFIXES include
)
find_library(BAT_LIBRARY
  NAMES BAT
  PATHS ${BAT_ROOT} $ENV{BATINSTALLDIR}
  PATH_SUFFIXES lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BAT
  FOUND_VAR BAT_FOUND
  REQUIRED_VARS
    BAT_LIBRARY
    BAT_INCLUDE_DIR
  VERSION_VAR BAT_VERSION
)
