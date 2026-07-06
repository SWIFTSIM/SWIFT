#!/bin/bash

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
#  Basic checks without parallel-HDF5.
#
#  Just the way we do this. No special significance. For sourcing from a check
#  script with modules loaded and the setup.sh script also sourced. Assumes
#  clean sources that have been had ./autogen.sh ran.
#
#  Peter W. Draper 20-JAN-2026.
#-
echo
echo "-------------"
echo "non-MPI build"
echo "-------------"
do_configure --disable-optimization --enable-debug
do_make
do_make check VERBOSE=1
do_make clean

echo
echo "---------"
echo "MPI build"
echo "---------"
do_configure --with-parmetis --disable-optimization
do_make
do_make check VERBOSE=1
do_make clean

echo
echo "----------------"
echo "No vectorisation"
echo "----------------"
do_configure --disable-vec
do_make
do_make check VERBOSE=1
do_make clean

echo
echo "-------------------"
echo "Basic Vectorisation"
echo "-------------------"
do_configure
do_make
do_make check VERBOSE=1
do_make clean

