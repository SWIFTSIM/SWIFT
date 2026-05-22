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
source ci/COSMA/intel-modules.sh
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

#  Formatting only from the tests.
do_make check TESTS=testFormat.sh
do_make clean

# exit

echo
echo "-----------"
echo "EAGLE build"
echo "-----------"
do_configure --with-subgrid=EAGLE --with-hydro=sphenix --disable-hand-vec
do_make
do_make clean

echo
echo "-----------"
echo "EAGLE build"
echo "-----------"
do_configure --with-subgrid=SPIN_JET_EAGLE --with-hydro=sphenix --disable-hand-vec
do_make
do_make clean

echo
echo "--------------"
echo "FLAMINGO build"
echo "--------------"
do_configure --with-subgrid=FLAMINGO --with-hydro=sphenix --disable-hand-vec
do_make
do_make clean

echo
echo "----------------------"
echo "GEAR + grackle 0 build"
echo "----------------------"
# Load grackle
module load grackle-swift/3.3.dev1

do_configure --with-subgrid=GEAR --with-hydro=sphenix --disable-hand-vec --with-grackle=${GRACKLE_HOME}/lib
do_make
do_make clean

echo
echo "----------------------"
echo "GEAR + grackle 3 build"
echo "----------------------"
do_configure --with-subgrid=GEAR-G3 --with-hydro=sphenix --disable-hand-vec --with-grackle=${GRACKLE_HOME}/lib
do_make
do_make clean

echo
echo "----------------"
echo "Prep1 loop build"
echo "----------------"
do_configure --with-hydro=sphenix --disable-hand-vec --with-stars=basic CFLAGS="-DEXTRA_STAR_LOOPS_1"
do_make
do_make clean

echo
echo "-----------------------------------------------"
echo "Prep1 + 2 loop build (similar to EAGLE-kinetic)"
echo "-----------------------------------------------"
do_configure --with-hydro=sphenix --disable-hand-vec --with-stars=basic CFLAGS="'-DEXTRA_STAR_LOOPS_1 -DEXTRA_STAR_LOOPS_2'"
do_make
do_make clean

echo
echo "--------------------------------------------------------"
echo "Prep2 + 3 loop build (similar to GEAR-mechanical mode 1)"
echo "--------------------------------------------------------"
do_configure --with-hydro=sphenix --disable-hand-vec --with-stars=basic CFLAGS="'-DEXTRA_STAR_LOOPS_2 -DEXTRA_STAR_LOOPS_3"
do_make
do_make clean

echo
echo "------------------------------------------------------------"
echo "Prep2 + 3 + 4 loop build (similar to GEAR-mechanical mode 2)"
echo "------------------------------------------------------------"
do_configure --with-hydro=sphenix --disable-hand-vec --with-stars=basic CFLAGS="'-DEXTRA_STAR_LOOPS_2 -DEXTRA_STAR_LOOPS_3 -DEXTRA_STAR_LOOPS_4'"
do_make
do_make clean

echo
echo "----------------------------"
echo "Prep1 + 2 + 3 + 4 loop build"
echo "----------------------------"
do_configure --with-hydro=sphenix --disable-hand-vec --with-stars=basic CFLAGS="-DEXTRA_STAR_LOOPS_1 -DEXTRA_STAR_LOOPS_2 -DEXTRA_STAR_LOOPS_3 -DEXTRA_STAR_LOOPS_4'"
do_make
do_make clean

#  Keep simple, may have a number of these happening.
exit
