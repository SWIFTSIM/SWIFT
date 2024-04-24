#!/bin/bash
set -o xtrace

# Resolution
N_label=n60

# Copy or download the initial conditions if they are not present
if [ ! -e demo_impact_"$N_label".hdf5 ]
then
    if [ -e ../DemoImpactInitCond/demo_impact_"$N_label".hdf5 ]
    then
        echo "Copying initial conditions from the DemoInitCond example..."
        cp ../DemoImpactInitCond/demo_impact_"$N_label".hdf5 ./
    else
        echo "Downloading initial conditions..."
        ./get_init_cond.sh
    fi
fi

# Download equation of state tables if not already present
if [ ! -e ../EoSTables/ANEOS_forsterite_S19.txt ]
then
    cd ../EoSTables
    ./get_eos_tables.sh
    cd -
fi

# Run SWIFT
../../../swift --hydro --self-gravity --threads=28 demo_impact_"$N_label".yml 2>&1 | tee output_"$N_label".txt

# Plot the snapshots
python3 plot_snapshots.py
