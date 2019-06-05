#!/bin/bash

# This example is based on the AGORA disk article (DOI: 10.3847/1538-4357/833/2/202)

# currently only the low resolution is available
sim=low

# make run.sh fail if a subcommand fails
set -e

# Generate the initial conditions if they are not present.
if [ ! -e $sim.hdf5 ]
then
    echo "Fetching initial glass file for the Sedov blast example..."
    ./getIC.sh $sim
fi

# Get the Grackle cooling table
if [ ! -e CloudyData_UVB=HM2012.h5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    ../../Cooling/getGrackleCoolingTable.sh 
fi

# copy the initial conditions
cp $sim.hdf5 agora_disk.hdf5
# Update the particle types
python3 changeType.py agora_disk.hdf5

# Run SWIFT
../../swift --cooling --hydro --self-gravity --threads=4 agora_disk.yml 2>&1 | tee output.log


echo "Changing smoothing length to be Gadget compatible"
python3 cleanupSwift.py agora_disk_0000.hdf5 agora_disk_IC.hdf5
python3 cleanupSwift.py agora_disk_0050.hdf5 agora_disk_500Myr.hdf5

echo "Fetching GEAR solution..."
./getSolution.sh


echo "Plotting..."
python3 plotSolution.py
