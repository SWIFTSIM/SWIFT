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
#  Checks that the "make dist" target works and we can build out of tree
#  as all good autotools projects should.
#
#  For sourcing from a check script with modules loaded and the setup.sh
#  script also sourced. Assumes clean sources that have been had ./autogen.sh
#  ran.
#
#  Peter W. Draper 20-JAN-2026.
#-
echo
echo "------------------------"
echo "Testing make dist target"
echo "------------------------"
do_configure --with-parmetis --disable-optimization
do_make dist
do_make distclean

echo "## Building from distribution tarball"
mkdir build
cd build
tar zxvf ../*.tar.gz
cd swift-*
do_configure --with-parmetis --disable-optimization
do_make
do_make check  VERBOSE=1
do_make clean
cd ../../
git clean -fdx

echo "-----------------"
echo "Out of tree build"
echo "-----------------"
./autogen.sh
mkdir build
cd build
../configure --with-parmetis --disable-optimization
do_make
do_make check VERBOSE=1
do_make clean
cd ../
git clean -fdx
./autogen.sh
