#!/bin/bash

#  Checks that the "make dist" target works and we can build out of tree
#  as all good autotools projects should.
#
#  For sourcing from a check script with modules loaded and the setup.sh
#  script also sourced. Assumes clean sources that have been had ./autogen.sh
#  ran.
#
#  Peter W. Draper 20-JAN-2026.

echo
echo "------------------------"
echo "Testing make dist target"
echo "------------------------"
do_make dist
do_make distclean

echo "Building from distribution tarball"
mkdir build
cd build
tar zxvf ../*.tar.gz
cd swift-*
../configure --with-parmetis --disable-optimization
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

