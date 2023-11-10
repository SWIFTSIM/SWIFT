#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e FastRotor_LR.hdf5 ]
then
    echo "Fetching Glass Files..."
    ./getGlass.sh
    echo "Generating the ICs"
    python ./makeIC_VP.py 
fi

# Run SWIFT
../../../swift --hydro --threads=16 ../FastRotor.yml 2>&1 > out.log 

# Plot the temperature evolution
python3 ./plot_all.py 0 60
