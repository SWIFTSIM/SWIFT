#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e glassPlane_128.hdf5 ]
then
    echo "Fetching initial glass file for the Yee vortex example..."
    ./getGlass.sh
fi
if [ ! -e yeeVortex.hdf5 ]
then
    echo "Generating initial conditions for the Yee vortex example..."
    python3 makeIC.py 128
fi

# Run SWIFT (configure with --with-hydro-dimension=2 --with-adiabatic-index=7/5)
../../../swift --hydro --threads=4 yeeVortex.yml 2>&1 | tee output.log

# Plot the solution at t = 10 against the analytic (stationary) vortex
python3 plotSolution.py 4
