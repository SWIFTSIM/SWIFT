#!/bin/bash -l

# This file is part of SWIFT.
# Copyright (C) 2026 p.w.draper@durham.ac.uk.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#+
#  Locally provided toolchain compilation and testing.
#
#  Peter W. Draper 20-JAN-2026.
#-
source ci/setup.sh

#  Start from clean sources.
git clean -fdx

#  And off we go.
./autogen.sh

#  Basic checks.
source ci/swift-checks.sh

#  More extensive tests that will use parallel if HDF5 available.
source ci/swift-checks-parallel.sh

#  Checks that distribution still works.
source ci/swift-checks-distribution.sh

exit
