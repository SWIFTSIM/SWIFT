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
#  Tests to be ran when changes are committed to a merge request.
#
#  Peter W. Draper 20-JAN-2026.
#-
#  Build toolchain.
source ci/intel-modules.sh
source ci/setup.sh

#  Clean sources.
git clean -fdx

#  And off we go.
./autogen.sh

echo
echo "-------------"
echo "non-MPI build"
echo "-------------"
do_configure --disable-optimization
do_make
do_make clean

echo
echo "---------"
echo "MPI build"
echo "---------"
do_configure --with-parmetis --disable-optimization
do_make
do_make clean

#  Keep simple, may have a number of these happening.
exit
