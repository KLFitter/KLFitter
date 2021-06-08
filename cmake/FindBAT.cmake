# Copyright (c) 2009--2021, the KLFitter developer team
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
# This file assumes that a properly installed version of BAT comes with
# pkgconfig files. Thus, use the pkgconfig module to locate BAT. Otherwise, if
# BAT is installed within the KLFitter build environment, the $BAT_ROOT
# environment variable will be set and can be used as a hint to locate the BAT
# installation. This script will define the following variables:
#
# BAT_FOUND    - True if the system has the BAT library
#
# BAT_VERSION - The version of the BAT library which was found
#
# and the following imported targets:
#
#   BAT::BAT   - The BAT library

find_package(PkgConfig)
pkg_check_modules(PC_BAT QUIET bat)

find_path(
  BAT_INCLUDEDIR
  NAMES BAT/BCLog.h BAT/BCMath.h
  PATHS ${PC_BAT_INCLUDE_DIRS} ${BAT_ROOT}
  PATH_SUFFIXES include)
find_library(
  BAT_LIBRARY
  NAMES BAT BATmodels BATmtf BATmvc
  PATHS ${PC_BAT_LIBRARY_DIRS} ${BAT_ROOT}
  PATH_SUFFIXES lib)

set(BAT_VERSION ${PC_BAT_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  BAT
  FOUND_VAR BAT_FOUND
  REQUIRED_VARS BAT_LIBRARY BAT_INCLUDEDIR
  VERSION_VAR BAT_VERSION)
