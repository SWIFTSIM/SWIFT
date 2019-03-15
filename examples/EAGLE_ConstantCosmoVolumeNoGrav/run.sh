#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e gravity_glassCube_32.hdf5 ]
then
    echo "Fetching initial gravity glass file for the constant cosmological box example..."
    ./getGlass.sh
fi

# Generate the initial conditions if they are not present.
if [ ! -e coolingtables ]
then
    echo "Fetching EAGLE Cooling Tables"
    ./getCoolingTable.sh
fi
if [ ! -e constantBox.hdf5 ]
then
    echo "Generating initial conditions for the uniform cosmo box example..."
    python3 makeIC.py
fi

# Run SWIFT
../swift --hydro --cosmology --cooling --threads=8 eagle_cooling_box.yml 2>&1 | tee output.log

# Plot the result
python3 plot_thermal_history.py eagle_cooling_box thermal_history_default_gadget2.pdf
