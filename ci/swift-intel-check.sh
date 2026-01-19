#!/bin/bash -l

#  SWIFT Intel toolchain complete set of builds and tests.


# When exiting in error report current configuration options.
function ONEXIT {
   if test "$?" != 0; then
      echo "Current configuration: $(grep "\./configure" config.log)"
   fi
}
trap ONEXIT EXIT

#  Lots of output.
set -e
set -x

#  Build toolchain.
source ci/intel-modules.sh

#  Clean sources.
git clean -fdx

#  And off we go.
./autogen.sh

echo
echo "non-MPI build"
echo "-------------"
./configure --disable-optimization --enable-debug
make -j 2
make check VERBOSE=1
make clean

echo
echo "MPI build"
echo "---------"
./configure --with-parmetis --disable-optimization
make -j 2
make check VERBOSE=1
make clean

echo
echo "No vectorisation"
echo "----------------"
./configure --disable-vec
make -j 2
make check VERBOSE=1
make clean

echo
echo "Basic Vectorisation"
echo "-------------------"
./configure
make -j 2
make check VERBOSE=1
make clean

echo
echo "Parallel HDF5 build"
echo "-------------------"
source ci/intel-modules-parallel.sh
./configure --with-parmetis --disable-optimization
make -j 2
make check VERBOSE=1
make clean

echo
echo "Testing dist target"
echo "-------------------"
make dist
make distclean

#  Try building it.
mkdir build
cd build
tar zxvf ../*.tar.gz
cd swift-*
./configure --with-parmetis --disable-optimization
make -j 2
make check VERBOSE=1
make clean
cd ../../
git clean -fdx

# Check that the debugging routines compile
./autogen.sh
./configure --with-parmetis --disable-vec --disable-optimization --enable-debugging-checks
make -j 2
make clean

# Check that the task-debugging output compiles
./configure --with-parmetis --disable-vec --disable-optimization --enable-task-debugging
make -j 2 
make clean

#  Add simple checks of gravity schemes.

./configure --with-parmetis --with-gravity=basic
make -j 2
make clean

./configure --with-parmetis --disable-gravitational-potential
make -j 2
make clean

#  Add simple checks of hydro schemes.

# GIZMO-MFM - Exact Riemann solver
./configure --with-parmetis --disable-vec --disable-optimization --with-hydro=gizmo-mfm --with-riemann-solver=exact
make -j 2
make clean

# GIZMO-MFV - HLLC Riemann solver
./configure --with-parmetis --disable-vec --disable-optimization --with-hydro=gizmo-mfv --with-riemann-solver=hllc
make -j 2
make clean

# MINIMAL
./configure --with-parmetis --disable-vec --disable-optimization --with-hydro=minimal
make -j 2
make clean

# PRESSURE-ENERGY
./configure --with-parmetis --disable-vec --disable-optimization --with-hydro=pressure-energy
make -j 2
make clean

# PRESSURE-ENTROPY
./configure --with-parmetis --disable-vec --disable-optimization --with-hydro=pressure-entropy
make -j 2
make clean

# PLANETARY PHYSICS
./configure --with-hydro=planetary --with-equation-of-state=planetary --enable-debugging-checks
make -j 2
make clean

# External gravity.

# Pointmass.
./configure --with-parmetis --disable-vec --disable-optimization --with-ext-potential=point-mass
make -j 2
make clean

# Disc patch.
./configure --with-parmetis --disable-vec --disable-optimization --with-ext-potential=disc-patch
make -j 2
make clean

# Isothermal.
./configure --with-parmetis --disable-vec --disable-optimization --with-ext-potential=isothermal
make -j 2
make clean

# Cooling functions.

# Const du
./configure --with-parmetis --disable-vec --disable-optimization --with-cooling=const-du
make -j 2
make clean

# Const lambda
./configure --with-parmetis --disable-vec --disable-optimization --with-cooling=const-lambda
make -j 2
make clean

# Const lambda
./configure --with-parmetis --disable-vec --disable-optimization --with-cooling=EAGLE --with-chemistry=EAGLE
make -j 2
make clean

# Full EAGLE-like stack
./configure --with-subgrid=EAGLE --with-hydro=sphenix --disable-hand-vec
make -j 2
make clean

# Stand-alone FoF
./configure --enable-stand-alone-fof
make -j 2
make clean

# Full EAGLE-XL-like stack
./configure --with-subgrid=EAGLE-XL --with-hydro=sphenix --disable-hand-vec
make -j 2
make clean

# Full EAGLE-SPIN-JET stack
./configure --with-subgrid=SPIN_JET_EAGLE-XL --with-hydro=sphenix --disable-hand-vec
make -j 2
make clean

exit
