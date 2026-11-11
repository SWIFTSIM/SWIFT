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

# Run SWIFT (configure with --with-hydro-dimension=2 --with-adiabatic-index=7/5).
# SWIFT_ARGS passes extra arguments, e.g. MAGMA's recommended
# SWIFT_ARGS="-P SPH:resolution_eta:1.5 -P SPH:CFL_condition:0.2"
../../../swift --hydro --threads=4 $SWIFT_ARGS yeeVortex.yml 2>&1 | tee output.log

# Plot the solution at t = 10 against the analytic (stationary) vortex
python3 plotSolution.py 4
