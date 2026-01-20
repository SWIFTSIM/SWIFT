#!/bin/bash -l

#  GCC toolchain compilation and testing.
source ci/gcc-modules.sh
source ci/setup.sh

#  Start from clean sources.
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
source ci/gcc-modules-parallel.sh

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
do_make check  VERBOSE=1
do_make clean
cd ../../
git clean -fdx

#  Do an out of tree build.
./autogen.sh
mkdir build
cd build
../configure --with-parmetis --disable-optimization
do_make
do_make check VERBOSE=1
do_make clean
cd ../
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

# Check that the gravity force checks compiles
./configure  --enable-gravity-force-checks=1000
do_make
do_make clean

#  Add simple checks of hydro schemes.

# GIZMO exact solver
./configure --with-parmetis --disable-vec --disable-optimization --with-hydro=gizmo-mfm --with-riemann-solver=exact
do_make
do_make clean

# GIZMO hllc solver
./configure --with-parmetis --disable-vec --disable-optimization --with-hydro=gizmo-mfm --with-riemann-solver=hllc
do_make
do_make clean

# MINIMAL
./configure --with-parmetis --disable-vec --disable-optimization --with-hydro=minimal
do_make
do_make clean

# ANARCHY
./configure --with-hydro=phantom
do_make
do_make check VERBOSE=1
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

# EAGLE
./configure --with-parmetis --disable-vec --disable-optimization --with-cooling=EAGLE --with-chemistry=EAGLE
do_make
do_make clean

# Chemistry models.

# Gear
./configure --disable-vec --with-chemistry=GEAR_10
do_make
do_make clean

# EAGLE
./configure --disable-vec --with-chemistry=EAGLE
do_make
do_make clean

# Velociraptor compliation check.
./configure --enable-dummy-velociraptor
do_make
do_make clean

# Stand-alone FoF
./configure --enable-stand-alone-fof
do_make
do_make clean

# Full EAGLE-like stack
./configure --with-subgrid=EAGLE --with-hydro=sphenix --disable-hand-vec
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
