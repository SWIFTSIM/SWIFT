#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "You need to provide the resolution (e.g. ./run.sh low)."
    echo "The possible options are low, med and high."
    exit
fi


# Generate the initial conditions if they are not present.
if [ ! -e $1.hdf5 ]
then
    echo "Fetching initial glass file for the Sedov blast example..."
    ./getIC.sh $1.hdf5
fi

cp $1.hdf5 agora_disk.hdf5
./change_type.py agora_disk.hdf5

# Run SWIFT
../swift -s -G -t 8 agora_disk.yml 2>&1 | tee output.log

# Plot the solution
# python plotSolution.py 5
