#!/bin/bash -l

#  SWIFT Intel toolchain complete set of builds and tests.

#  Build toolchain.
source ci/intel-modules.sh
source ci/setup.sh

#  Clean sources.
git clean -fdx

#  And off we go.
./autogen.sh

echo
echo "non-MPI build"
echo "-------------"
./configure --disable-optimization --enable-debug
do_make
do_make check VERBOSE=1
do_make clean

echo
echo "MPI build"
echo "---------"
./configure --with-parmetis --disable-optimization
do_make
do_make check VERBOSE=1
do_make clean

echo
echo "No vectorisation"
echo "----------------"
./configure --disable-vec
do_make
do_make check VERBOSE=1
do_make clean

echo
echo "Basic Vectorisation"
echo "-------------------"
./configure
do_make
do_make check VERBOSE=1
do_make clean

echo
echo "Parallel HDF5 build"
echo "-------------------"
source ci/intel-modules-parallel.sh
./configure --with-parmetis --disable-optimization
do_make
do_make check VERBOSE=1
do_make clean

echo
echo "Testing dist target"
echo "-------------------"
do_make dist
do_make distclean

#  Try building it.
mkdir build
cd build
tar zxvf ../*.tar.gz
cd swift-*
./configure --with-parmetis --disable-optimization
do_make
do_make check VERBOSE=1
do_make clean
cd ../../
git clean -fdx

# Check that the debugging routines compile
./autogen.sh
./configure --with-parmetis --disable-vec --disable-optimization --enable-debugging-checks
do_make
do_make clean

# Check that the task-debugging output compiles
./configure --with-parmetis --disable-vec --disable-optimization --enable-task-debugging
do_make
do_make clean

#  Add simple checks of gravity schemes.

./configure --with-parmetis --with-gravity=basic
do_make
do_make clean

./configure --with-parmetis --disable-gravitational-potential
do_make
do_make clean

#  Add simple checks of hydro schemes.

# GIZMO-MFM - Exact Riemann solver
./configure --with-parmetis --disable-vec --disable-optimization --with-hydro=gizmo-mfm --with-riemann-solver=exact
do_make
do_make clean

# GIZMO-MFV - HLLC Riemann solver
./configure --with-parmetis --disable-vec --disable-optimization --with-hydro=gizmo-mfv --with-riemann-solver=hllc
do_make
do_make clean

# MINIMAL
./configure --with-parmetis --disable-vec --disable-optimization --with-hydro=minimal
do_make
do_make clean

# PRESSURE-ENERGY
./configure --with-parmetis --disable-vec --disable-optimization --with-hydro=pressure-energy
do_make
do_make clean

# PRESSURE-ENTROPY
./configure --with-parmetis --disable-vec --disable-optimization --with-hydro=pressure-entropy
do_make
do_make clean

# PLANETARY PHYSICS
./configure --with-hydro=planetary --with-equation-of-state=planetary --enable-debugging-checks
do_make
do_make clean

# External gravity.

# Pointmass.
./configure --with-parmetis --disable-vec --disable-optimization --with-ext-potential=point-mass
do_make
do_make clean

# Disc patch.
./configure --with-parmetis --disable-vec --disable-optimization --with-ext-potential=disc-patch
do_make
do_make clean

# Isothermal.
./configure --with-parmetis --disable-vec --disable-optimization --with-ext-potential=isothermal
do_make
do_make clean

# Cooling functions.

# Const du
./configure --with-parmetis --disable-vec --disable-optimization --with-cooling=const-du
do_make
do_make clean

# Const lambda
./configure --with-parmetis --disable-vec --disable-optimization --with-cooling=const-lambda
do_make
do_make clean

# Const lambda
./configure --with-parmetis --disable-vec --disable-optimization --with-cooling=EAGLE --with-chemistry=EAGLE
do_make
do_make clean

# Full EAGLE-like stack
./configure --with-subgrid=EAGLE --with-hydro=sphenix --disable-hand-vec
do_make
do_make clean

# Stand-alone FoF
./configure --enable-stand-alone-fof
do_make
do_make clean

# Full EAGLE-XL-like stack
./configure --with-subgrid=EAGLE-XL --with-hydro=sphenix --disable-hand-vec
do_make
do_make clean

# Full EAGLE-SPIN-JET stack
./configure --with-subgrid=SPIN_JET_EAGLE-XL --with-hydro=sphenix --disable-hand-vec
do_make
do_make clean

exit
