#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e nfw.hdf5 ]
then
    echo "Fetching initial conditions file for the example..."
    ./getICs.sh
fi


# Run SWIFT
../../../swift --hydro --external-gravity --self-gravity --threads=14  NFW_Hydrostatic.yml

# Compute gas density profile
python3  plotGasDensityProfile.py
