#!/bin/bash

# Get the initial conditions if they are not present.
if [ ! -e uranus_1e6.hdf5 ]
then
    echo "Fetching initial conditions file for the Uranus impact example..."
    ./get_init_cond.sh
fi

# Get the EoS tables if they are not present.
cd ../EoSTables
if [ ! -e HM80_HHe.txt ] || [ ! -e HM80_ice.txt ] || [ ! -e HM80_rock.txt ] 
then
    echo "Fetching equations of state tables for the Uranus impact example..."
    ./get_eos_tables.sh
fi
cd ../GiantImpact

# Run SWIFT
../../swift -s -G -t 8 uranus_1e6.yml 2>&1 | tee output.log

# Plot the solution
python3 plot_solution.py
