#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e glassCube_32.hdf5 ]
then
    echo "Fetching initial glass file for the Supernovae feedback example..."
    ./getGlass.sh
fi
if [ ! -e stellar_evolution.hdf5 ]
then
    echo "Generating initial conditions for the 3D stellar evolution example..."
    python makeIC.py
fi

# Get the Yield tables
if [ ! -e yieldtables ]
then
    echo "Fetching Yield tables..."
    ./getYieldTable.sh
fi


../../swift --limiter --feedback --stars --hydro --external-gravity --threads=4 stellar_evolution.yml 2>&1 | tee output.log

#python check_stellar_evolution.py
python plot_box_evolution.py
