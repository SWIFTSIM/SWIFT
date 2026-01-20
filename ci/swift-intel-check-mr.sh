#!/bin/bash -l

#  Tests to be ran when changes are committed to a merge request.

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
