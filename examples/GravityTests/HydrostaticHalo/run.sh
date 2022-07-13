#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e Hydrostatic.hdf5 ]
then
    echo "Generating initial conditions for the isothermal potential box example..."
    python3 makeIC.py 100000
fi

# Run for 10 dynamical times
../../../swift --external-gravity --hydro --threads=1 hydrostatic.yml 2>&1 | tee output.log

echo "Plotting density profiles"
mkdir plots
mkdir plots/density_profile
python3 density_profile.py 2. 200 300

echo "Plotting internal energy profiles"
mkdir plots/internal_energy
python3 internal_energy_profile.py 2. 200 300

echo "Plotting radial velocity profiles"
mkdir plots/radial_velocity_profile
python3 velocity_profile.py 2. 200 300

echo "Plotting energy as a function of time"
python3 test_energy_conservation.py 300
