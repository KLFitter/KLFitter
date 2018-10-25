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

# Helper macro for building the project's executables.
macro( KLFitter_add_executable name )
  add_executable( ${name} ${ARGN} )
  target_link_libraries( ${name} KLFitter )
  install( TARGETS ${name}
    EXPORT KLFitterTargets
    RUNTIME DESTINATION bin )
endmacro( KLFitter_add_executable )

# Build the example executable(s) of the project.
KLFitter_add_executable( example-top-ljets.exe example-top-ljets.cxx )
