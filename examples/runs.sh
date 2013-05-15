#!/bin/bash

# Set some global stuff
export OMP_WAIT_POLICY=PASSIVE

# Loop over number of cores
for cpu in {1..32}
do

    # Sod-Shock runs
    if [ ! -e SodShock_${cpu}.dump ]
    then
        ./test -c 1.0 -t $cpu -f sodShock.hdf5 -m 0.1 -w 5000 -d 1.0 > SodShock_${cpu}.dump
    fi
    if [ ! -e SodShock_fixed_${cpu}.dump ]
    then
        ./test_fixdt -r 1000 -t $cpu -f sodShock.hdf5 -m 0.1 -w 5000 -d 1e-3 > SodShock_fixed_${cpu}.dump
    fi
    
    # Sedov blast
    if [ ! -e SedovBlast_${cpu}.dump ]
    then
        ./test -c 0.2 -t $cpu -f sedov.hdf5 -m 1.0 -w 5000 -d 0.2 > SedovBlast_${cpu}.dump
    fi
    if [ ! -e SedovBlast_fixed_${cpu}.dump ]
    then
        ./test_fixdt -r 4096 -t $cpu -f sedov.hdf5 -m 1.0 -w 5000 -d 5e-5 > SedovBlast_fixed_${cpu}.dump
    fi
    
    # Cosmological volume
    if [ ! -e CosmoVolume_${cpu}.dump ]
    then
        ./test -c 0.01 -t $cpu -f cosmoVolume.hdf5 -m 0.6 -w 5000 -d 0.01 > CosmoVolume_${cpu}.dump
    fi
    if [ ! -e CosmoVolume_fixed_${cpu}.dump ]
    then
        ./test_fixdt -r 256 -t $cpu -f cosmoVolume.hdf5 -m 0.6 -w 5000 -d 1e-8 > CosmoVolume_fixed_${cpu}.dump
    fi

done

