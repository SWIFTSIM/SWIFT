#!/bin/bash -l

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
#  Runtime checks of SWIFT. So we build the binaries and do some small checks
#  that this works for some of the examples.
#
#  Note that we no longer use the -e flag to check floating point operations
#  as as that is no longer reliable with Intel 2024 onwards.
#
#  Peter W. Draper 20-JAN-2026.
#-
#  Build toolchain.
source ci/COSMA/intel-modules.sh
source ci/setup.sh

#  Extra modules for runtime.
module load grackle-swift/3.3.dev1

#  Need 2.70 for Fortran support.
module load autoconf

#  Clean sources.
git clean -fdx

#  And off we go.
./autogen.sh

echo 
echo "------------------------"
echo "Building SWIFT binaries."
echo "------------------------"
do_configure --with-parmetis --enable-debugging-checks --disable-vec --enable-debug --with-ext-potential=point-mass
do_make

#  Now run some simple checks.
cd examples/GravityTests/ExternalPointMass
echo "------------------------"
echo "External point mass test, 1 thread"
echo "------------------------"
python makeIC.py 10000
do_run ../../../swift -g -t 1 externalPointMass.yml

echo "---------------------------------------------"
echo "External point mass test, 1 thread, drift all"
echo "---------------------------------------------"
do_run ../../../swift -g -D -t 1 externalPointMass.yml

#  Test all particle types with the EAGLE-6 model
cd ../../EAGLE_low_z/EAGLE_6
echo "-----------------------"
echo "EAGLE-6 test, 4 threads"
echo "-----------------------"

#  Avoid downloading it is not necessary on COSMA which has direct access.
#wget http://virgodb/swift-webstorage/ICs/EAGLE_low_z/EAGLE_ICs_6.hdf5
ln -s /cosma5/data/Swift/web-storage/ICs/EAGLE_low_z/EAGLE_ICs_6.hdf5 EAGLE_ICs_6.hdf5
do_run ../../../swift -c -s -S -G -t 4 eagle_6.yml -n 16

# Check also a dry-run
echo "--------------------------------"
echo "EAGLE-6 test, 4 threads, dry-run"
echo "--------------------------------"
do_run ../../../swift -c -s -S -G -t 4 eagle_6.yml -n 16 -d

# Proper multi-threaded hydro-test.
cd ../../HydroTests/SodShock_3D
echo "--------------------------------------------"
echo "SodShock threaded test, 4 threads, 256 steps"
echo "--------------------------------------------"
#wget http://virgodb/swift-webstorage/ICs/glassCube_64.hdf5
#wget http://virgodb/swift-webstorage/ICs/glassCube_32.hdf5
ln -s /cosma5/data/Swift/web-storage/ICs/glassCube_64.hdf5 glassCube_64.hdf5
ln -s /cosma5/data/Swift/web-storage/ICs/glassCube_32.hdf5 glassCube_32.hdf5
python makeIC.py
do_run ../../../swift -s -t 4 -n 256 sodShock.yml

#  Bigger test, requires a lot of resources, so don't run as long.
unset I_MPI_HYDRA_BOOTSTRAP
echo "------------------------------------------------"
echo "SodShock MPI test, 16 ranks, 4 threads, 64 steps"
echo "------------------------------------------------"
do_run mpirun -np 16 ../../../swift_mpi -s -t 4 -n 64 sodShock.yml -PScheduler:max_top_level_cells:24
cd ../../../

echo 
echo "-----------------------"
echo "Building SWIFT binaries"
echo "-----------------------"
do_make clean
do_configure --with-parmetis --enable-debugging-checks --enable-debug
do_make

# Test EAGLE 12 with external gravity over MPI.
cd examples/EAGLE_low_z/EAGLE_12
echo "------------------------------------------------------------------"
echo "EAGLE_12 MPI test with ext. gravity, 4 ranks, 8 threads, 128 steps"
echo "------------------------------------------------------------------"
#wget http://virgodb/swift-webstorage/ICs/EAGLE_low_z/EAGLE_ICs_12.hdf5
ln -s /cosma5/data/Swift/web-storage/ICs/EAGLE_low_z/EAGLE_ICs_12.hdf5 EAGLE_ICs_12.hdf5
do_run mpirun -np 4 ../../../swift_mpi -g -s -t 8 -n 1024 -PRestarts:onexit:1 eagle_12.yml

#  Restart this.
echo "----------"
echo "Restarting"
echo "----------"
do_run mpirun -np 4 ../../../swift_mpi -g -s -t 8 -n 1200 eagle_12.yml -r
cd ../../../

echo "---------------------"
echo "1D check of SodShock."
echo "---------------------"
do_make clean
do_configure --with-hydro-dimension=1 --with-parmetis --enable-debugging-checks --disable-vec --enable-debug
do_make
cd examples/HydroTests/SodShock_1D
python makeIC.py
do_run ../../../swift -s -t 1 sodShock.yml
cd ../../../

echo "--------------"
echo "Sink particles"
do_make clean
echo "--------------"
echo "debugging-checks are disabled"
#do_configure --disable-mpi --with-chemistry=GEAR_10 --with-cooling=grackle_0 --with-stars=GEAR --with-star-formation=GEAR --with-feedback=GEAR --with-sink=GEAR --with-kernel=wendland-C2 --enable-debugging-checks --with-grackle=${GRACKLE_HOME}/lib
do_configure --disable-mpi --with-chemistry=GEAR_10 --with-cooling=grackle_0 --with-stars=GEAR --with-star-formation=GEAR --with-feedback=GEAR --with-sink=GEAR --with-kernel=wendland-C2 --with-grackle=${GRACKLE_HOME}/lib
do_make
cd examples/SinkParticles/PlummerSphere
ln -s /cosma5/data/Swift/web-storage/ICs/test_sink.hdf5
ln -s /cosma5/data/Swift/web-storage/CoolingTables/CloudyData_UVB=HM2012.h5
ln -s /cosma5/data/Swift/web-storage/FeedbackTables/POPIIsw.h5
do_run ../../../swift --hydro --sinks --stars --self-gravity --feedback --cooling --threads=4 params.yml

exit
