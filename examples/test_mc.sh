#!/bin/bash

# Set the default OpenMP behaviour
export OMP_WAIT_POLICY=PASSIVE

# Re-build the code on the target architecture
# cd ..
# ./configure --with-gcc-arch=native
# cd src
# make clean && make -j 8
# cd ../examples
# make clean && make -j 8

# Clear up the memory first
# ./memhog `free -g | grep Mem | awk '{print int(0.9*$2)}'`

# loop over number of CPUs
for cpu in {1..32}
do

    # Set some environment variables for OpenMP
    export OMP_NUM_THREADS=$cpu
    export OMP_THREAD_LIMIT=$cpu
    export OMP_PROC_BIND=TRUE
    
    # ./test -r 1000 -t $cpu -b "1400 1400 1400" -N 74240 -c small/Coordinates.txt.gz -s "700 700 700" -p 0 -h small/SmoothingLength.txt.gz -m 470 -z 200 > small_${cpu}.dump

    # ./test -r 100 -t $cpu -b "100 100 100" -N 3558892 -c snap_C09/Coordinates.txt.gz -s "50 50 50" -p 0 -h snap_C09/SmoothingLength.txt.gz -m 6.138 > snap_C09_${cpu}.dump

    ./test -r 100 -t $cpu -f scaling/snap_023_z000p503.hdf5 -m 0.5 -z 400 -w 5000 -d 1.0 > scaling_${cpu}.dump

    ./test -r 2000 -t $cpu -f SedovBlast/sedov.hdf5 -m 5e-1 -d 5.0e-4 > sedov_${cpu}.dump
    
    ./test -r 2000 -t $cpu -f PertubedBox/perturbedBox.hdf5 -m 0.05 -d 4e-3 > perturbed_${cpu}.dump

done
