#!/bin/bash

# currently only the low resolution is available
sim=low

# enable cooling or not
cooling=0

set -e

flag=

if [ $cooling -eq 1 ]
then
    flag=-C
fi

# Generate the initial conditions if they are not present.
if [ ! -e $sim.hdf5 ]
then
    echo "Fetching initial glass file for the Sedov blast example..."
    ./getIC.sh $sim.hdf5
fi

# Get the Grackle cooling table
if [ ! -e CloudyData_UVB=HM2012.h5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    ../getCoolingTable.sh
fi

cp $sim.hdf5 agora_disk.hdf5
python3 change_type.py agora_disk.hdf5

# Run SWIFT
#../swift $flag -s -G -t 4 agora_disk.yml 2>&1 | tee output.log


echo "Changing smoothing length to be Gadget compatible"
python cleanup_swift.py agora_disk_0000.hdf5 agora_disk_IC.hdf5
python cleanup_swift.py agora_disk_0050.hdf5 agora_disk_500Myr.hdf5

echo "Fetching Gear solution..."
./getSolution.sh $flag

echo "Plotting..."
python plot_solution.py $flag
