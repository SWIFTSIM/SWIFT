#!/bin/bash

set -eo pipefail

# Generate the initial conditions if they are not present.
if [ ! -e glassPlane_128.hdf5 ]
then
    echo "Fetching initial glass file for the Taylor-Green vortex example..."
    ./getGlass.sh
fi
if [ ! -e taylorGreen.hdf5 ]
then
    echo "Generating initial conditions for the Taylor-Green vortex example..."
    python3 makeIC.py
fi

# Run SWIFT (configure with --with-hydro-dimension=2)
# SWIFT_ARGS passes extra arguments, e.g. MAGMA's recommended
# SWIFT_ARGS="-P SPH:resolution_eta:1.5 -P SPH:CFL_condition:0.2"
../../../swift --hydro --threads=4 $SWIFT_ARGS taylorGreen.yml 2>&1 | tee output.log

# Plot the decay of the kinetic energy
python3 plotSolution.py
