#!/bin/bash -l

# This file is part of SWIFT.
# Copyright (C) 2026 schaller@strw.leidenuniv.nl
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
#  Runtime checks of SWIFT. So we build the binaries and do some small checks
#  that this works for some of the examples.
#
#  Note that we no longer use the -e flag to check floating point operations
#  as as that is no longer reliable with Intel 2024 onwards.
#
#  Matthieu Schaller 12-MAR-2026.
#-
#  Build toolchain.
source ci/COSMA/intel-modules.sh
source ci/setup.sh

#  Need 2.70 for Fortran support.
module load autoconf

#  Clean sources.
git clean -fdx

#  And off we go.
./autogen.sh

echo "--------------"
echo "SIDM          "
echo "--------------"
do_configure --with-sidm=Basic
do_make -j 2
do_make check
do_make clean

echo "--------------"
echo "SIDM + debugs "
echo "--------------"
do_configure --with-sidm=Basic --enable-debugging-checks
do_make -j 2
do_make check
do_make clean


echo "--------------"
echo "SIDM + EAGLE  "
echo "--------------"
do_configure --with-sidm=Basic --with-subgrid=EAGLE
do_make -j 2
do_make check
do_make clean


exit

