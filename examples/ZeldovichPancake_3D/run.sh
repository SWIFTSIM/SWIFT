#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e zeldovichPancake.hdf5 ]
then
    echo "Generating initial conditions for the 3D Zeldovich pancake example..."
    python makeIC.py
fi

# Run SWIFT
../swift -a -s -c -G -t 8 zeldovichPancake.yml 2>&1 | tee output.log

# Plot the result
for i in {0..119}
do 
    python plotSolution.py $i
done
