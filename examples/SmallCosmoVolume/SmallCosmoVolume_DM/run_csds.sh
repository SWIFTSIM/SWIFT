#!/bin/bash

#Path to velociraptor executable
velociraptor_path=/path/to/velociraptor

 # Generate the initial conditions if they are not present.
if [ ! -e small_cosmo_volume.hdf5 ]
then
    echo "Fetching initial conditions for the small cosmological volume example..."
    ./getIC.sh
fi

# Run SWIFT
../../../swift --cosmology --self-gravity --csds --threads=12 small_cosmo_volume_dm.yml 2>&1 | tee output.log

if [ ! -d "stf_output" ]
then
    mkdir stf_output
fi

# Run velociraptor
$velociraptor_path/stf -I 2 -C vrconfig_3dfof_subhalos_SO_hydro.cfg -i snap_0031 -o stf_output/snap_0031

# Make a movie of the evolution of the DM particles
python3 csds_analysis.py csds_index_0000.dump --halo 1
