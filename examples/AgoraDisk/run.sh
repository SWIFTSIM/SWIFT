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

# Get the Grackle cooling table
if [ ! -e CloudyData_UVB=HM2012.h5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    ../getCoolingTable.sh
fi

cp $1.hdf5 agora_disk.hdf5
./change_type.py agora_disk.hdf5

# Run SWIFT
../swift -C -s -G -t 4 agora_disk.yml 2>&1 | tee output.log
exit

echo "Changing smoothing length to be Gadget compatible"
./cleanup_swift.py agora_disk_0000.hdf5 agora_disk_IC.hdf5
./cleanup_swift.py agora_disk_0050.hdf5 agora_disk_500Myr.hdf5

if [ ! -e snapshot_0000 ] || [ ! -e snapshot_0500 ]
then
    echo "Fetching Gear solution..."
    ./getSolution.sh
fi

echo "Plotting..."
./plot_solution.py
