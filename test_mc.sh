#!/bin/bash

# Set the default OpenMP behaviour
export OMP_WAIT_POLICY=PASSIVE

# Clear up the memory first
# ./memhog `free -g | grep Mem | awk '{print int(0.9*$2)}'`

# loop over number of CPUs
for cpu in {1..32}
do

    # Set some environment variables for OpenMP
    export OMP_NUM_THREADS=$cpu
    export OMP_THREAD_LIMIT=$cpu
    export OMP_PROC_BIND=TRUE
    
    # ./test -r 100 -t $cpu -q $cpu -b "25000 25000 25000" -N 149646 -c data/Coordinates.txt -s "12500 12500 12500" -h data/SmoothingLength.txt > test_${cpu}.dump

    ./test -r 100 -t $cpu -q $cpu -b "100 100 100" -N 3664514 -c data2/Coordinates.txt -h data2/SmoothingLength.txt > test2_${cpu}.dump

    # ./test -r 100 -t $cpu -q $cpu -b "140 140 140" -N 7741 -c shrink/Coordinates.txt -s "70 70 70" -h shrink/SmoothingLength.txt > shrink_${cpu}.dump

done
