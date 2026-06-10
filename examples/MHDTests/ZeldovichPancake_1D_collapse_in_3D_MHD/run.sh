#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e zeldovichPancake.hdf5 ]
then
    echo "Generating initial conditions for the 1D Zeldovich pancake collapse example with magnetic field..."
    python3 makeIC_perpB.py
fi

# Run SWIFT
../../../swift --hydro --cosmology --self-gravity --threads=8 zeldovichPancake_MHD.yml 2>&1 | tee output.log

# Plot the result
for i in {0..119}
do 
    python3 plot_solution_perpB.py $i
done
