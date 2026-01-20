#!/bin/bash

#  Basic checks without parallel-HDF5.
#
#  Just the way we do this. No special significance. For sourcing from a check
#  script with modules loaded and the setup.sh script also sourced. Assumes
#  clean sources that have been had ./autogen.sh ran.
#
#  Peter W. Draper 20-JAN-2026.

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

