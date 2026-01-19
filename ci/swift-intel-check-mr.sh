#!/bin/bash -l

#  Tests to be ran when changes are committed to a merge request.

# When exiting in error report current configuration options.
function ONEXIT {
   if test "$?" != 0; then
      echo "Current configuration: $(grep "\./configure" config.log)"
   fi
}
trap ONEXIT EXIT

source intel-modules.sh
#module load intel_comp/2024.2.0 compiler-rt tbb compiler mpi
#module load hdf5/1.14.4
#module load fftw/3.3.10
#module load gsl/2.8
#module load parmetis/4.0.3
#module load python/3.12.4
#module load llvm/20.1.0
#export CLANG_FORMAT_CMD="clang-format"

set -e
set -x

git clean -fdx

./autogen.sh

echo
echo "non-MPI build"
echo "-------------"
./configure --disable-optimization
make
make clean

echo
echo "MPI build"
echo "---------"
./configure --with-parmetis --disable-optimization
make
make clean

#  Keep simple, may have a number of these happening.
exit


echo
echo "No vectorisation"
echo "----------------"
./configure --disable-vec
make
make clean

echo
echo "Basic Vectorisation"
echo "-------------------"
./configure
make
make clean

echo
echo "Logger"
echo "-------------------"
./configure --enable-logger
make
make clean

echo
echo "Low-memory model"
echo "-------------------"
./configure --with-hydro=none --with-stars=none --disable-gravitational-potential
make
make clean


echo
echo "Parallel HDF5 build"
echo "-------------------"
source intel-modules-parallel.sh
#module unload hdf5
#module load parallel_hdf5/1.14.4
./configure --with-parmetis --disable-optimization
make
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
make
make check VERBOSE=1
make clean
cd ../../
git clean -fdx

# Check that the debugging routines compile
./autogen.sh
./configure --with-parmetis --disable-vec --disable-optimization --enable-debugging-checks
make
make clean

# Check that the task-debugging output compiles
./configure --with-parmetis --disable-vec --disable-optimization --enable-task-debugging
make
make clean

#  Add simple checks of hydro schemes.

# GIZMO
./configure --with-parmetis --disable-vec --disable-optimization --with-hydro=gizmo --with-riemann-solver=exact
make
make clean

# MINIMAL
./configure --with-parmetis --disable-vec --disable-optimization --with-hydro=minimal
make
make clean

# DEFAULT
./configure --with-parmetis --disable-vec --disable-optimization --with-hydro=default
make
make clean

# PRESSURE
./configure --with-parmetis --disable-vec --disable-optimization --with-hydro=hopkins
make
make clean

# External gravity.

# Pointmass.
./configure --with-parmetis --disable-vec --disable-optimization --with-ext-potential=point-mass
make
make clean

# Disc patch.
./configure --with-parmetis --disable-vec --disable-optimization --with-ext-potential=disc-patch
make
make clean

# Isothermal.
./configure --with-parmetis --disable-vec --disable-optimization --with-ext-potential=isothermal
make
make clean

# Cooling functions.

# Const du
./configure --with-parmetis --disable-vec --disable-optimization --with-cooling=const-du
make
make clean

# Const lambda
./configure --with-parmetis --disable-vec --disable-optimization --with-cooling=const-lambda
make
make clean

# Full physics

# EAGLE
./configure --with-parmetis --disable-optimization --with-subgrid=EAGLE
make
make clean

# EAGLE-XL
./configure --with-parmetis --disable-optimization --with-subgrid=EAGLE-XL
make
make clean

# QLA
./configure --with-parmetis --disable-optimization --with-subgrid=QLA
make
make clean

exit
