#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e EAGLE_ICs_6.hdf5 ]
then
    echo "Fetching initial conditions for the EAGLE 6Mpc example..."
    ./getIC.sh
fi

# Adapt the initial conditions
python adaptIC.py

# Make first run for a hydrodynamical box
./configure --disable-hand-vec --with-hydro=sphenix --with-kernel=cubic-spline

./swift --hydro --threads=4 -n 16 density_test_Hydro.yml

# Make second run for a DMONLY box
./configure --with-sidm=yes

./swift --sidm -G --threads=4 -n 16 density_test_DM.yml

# Plot solution
python plot_solution.py
