#!/bin/bash

#  Basic checks with parallel-HDF5.
#
#  Just the way we do this. No special significance. For sourcing from a check
#  script with modules loaded and the setup.sh script also sourced. Assumes
#  clean sources that have been had ./autogen.sh ran.

echo
echo "-------------------"
echo "Parallel HDF5 build"
echo "-------------------"
do_configure --with-parmetis --disable-optimization
do_make
do_make check VERBOSE=1
do_make clean

echo "--------------------------------"
echo "Check debugging routines compile"
echo "--------------------------------"
do_configure --with-parmetis --disable-vec --disable-optimization --enable-debugging-checks
do_make
do_make clean

echo "------------------------------------"
echo "Check task-debugging output compiles"
echo "------------------------------------"
do_configure --with-parmetis --disable-vec --disable-optimization --enable-task-debugging
do_make
do_make clean

echo "-----------------------------------"
echo "Check gravity force checks compiles"
echo "-----------------------------------"
do_configure  --enable-gravity-force-checks=1000
do_make
do_make clean

echo "-------------------------------"
echo "simple checks of hydro schemes."
echo "-------------------------------"
echo "GIZMO exact solver"
do_configure --with-parmetis --disable-vec --disable-optimization --with-hydro=gizmo-mfm --with-riemann-solver=exact
do_make
do_make clean

echo "GIZMO hllc solver"
do_configure --with-parmetis --disable-vec --disable-optimization --with-hydro=gizmo-mfm --with-riemann-solver=hllc
do_make
do_make clean

echo "GIZMO-MFV - HLLC Riemann solver"
do_configure --with-parmetis --disable-vec --disable-optimization --with-hydro=gizmo-mfv --with-riemann-solver=hllc
do_make
do_make clean

echo "MINIMAL"
do_configure --with-parmetis --disable-vec --disable-optimization --with-hydro=minimal
do_make
do_make clean

echo "ANARCHY"
do_configure --with-hydro=phantom
do_make
do_make check VERBOSE=1
do_make clean

echo "PRESSURE-ENERGY"
./configure --with-parmetis --disable-vec --disable-optimization --with-hydro=pressure-energy
do_make
do_make clean

echo "PRESSURE-ENTROPY"
./configure --with-parmetis --disable-vec --disable-optimization --with-hydro=pressure-entropy
do_make
do_make clean

echo "PLANETARY PHYSICS"
./configure --with-hydro=planetary --with-equation-of-state=planetary --enable-debugging-checks
do_make
do_make clean

echo "----------------"
echo "External gravity"
echo "----------------"
echo "Pointmass"
do_configure --with-parmetis --disable-vec --disable-optimization --with-ext-potential=point-mass
do_make
do_make clean

echo "Disc patch"
do_configure --with-parmetis --disable-vec --disable-optimization --with-ext-potential=disc-patch
do_make
do_make clean

echo "Isothermal"
do_configure --with-parmetis --disable-vec --disable-optimization --with-ext-potential=isothermal
do_make
do_make clean

echo "-----------------"
echo "Cooling functions"
echo "-----------------"
echo "Const du"
do_configure --with-parmetis --disable-vec --disable-optimization --with-cooling=const-du
do_make
do_make clean

echo "Const lambda"
do_configure --with-parmetis --disable-vec --disable-optimization --with-cooling=const-lambda
do_make
do_make clean

echo "EAGLE"
do_configure --with-parmetis --disable-vec --disable-optimization --with-cooling=EAGLE --with-chemistry=EAGLE
do_make
do_make clean

echo "----------------"
echo "Chemistry models"
echo "----------------"
echo "Gear"
do_configure --disable-vec --with-chemistry=GEAR_10
do_make
do_make clean

echo "EAGLE"
do_configure --disable-vec --with-chemistry=EAGLE
do_make
do_make clean

echo "------------------------------"
echo "Velociraptor compilation check"
echo "------------------------------"
do_configure --enable-dummy-velociraptor
do_make
do_make clean

echo "---------------"
echo "Stand-alone FoF"
echo "---------------"
do_configure --enable-stand-alone-fof
do_make
do_make clean

echo "---------------------"
echo "Full EAGLE-like stack"
echo "---------------------"
do_configure --with-subgrid=EAGLE --with-hydro=sphenix --disable-hand-vec
do_make
do_make clean

echo "------------------------"
echo "Full EAGLE-XL-like stack"
echo "------------------------"
do_configure --with-subgrid=EAGLE-XL --with-hydro=sphenix --disable-hand-vec
do_make
do_make clean

echo "-------------------------"
echo "Full EAGLE-SPIN-JET stack"
echo "-------------------------"
do_configure --with-subgrid=SPIN_JET_EAGLE-XL --with-hydro=sphenix --disable-hand-vec
do_make
do_make clean
