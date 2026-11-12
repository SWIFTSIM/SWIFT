#!/bin/bash

set -eo pipefail

# Generate the initial conditions if they are not present.
if [ ! -e glassCube_32.hdf5 ]
then
    echo "Fetching initial glass file for the strong blast wave example..."
    ./getGlass.sh
fi
if [ ! -e blastWave.hdf5 ]
then
    echo "Generating initial conditions for the strong blast wave example..."
    python3 makeIC.py
fi

# Run SWIFT
# SWIFT_ARGS passes extra arguments, e.g. MAGMA's recommended
# SWIFT_ARGS="-P SPH:resolution_eta:1.5 -P SPH:CFL_condition:0.2"
../../../swift --hydro --threads=4 $SWIFT_ARGS blastWave.yml 2>&1 | tee output.log

# Plot the solution
python3 plotSolution.py 2
