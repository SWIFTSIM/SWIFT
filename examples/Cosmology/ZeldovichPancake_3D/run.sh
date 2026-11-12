#!/bin/bash

set -eo pipefail

# Generate the initial conditions if they are not present.
if [ ! -e zeldovichPancake.hdf5 ]
then
    echo "Generating initial conditions for the 3D Zeldovich pancake example..."
    python3 makeIC.py
fi

# Run SWIFT
# SWIFT_ARGS passes extra arguments, e.g. MAGMA's recommended
# SWIFT_ARGS="-P SPH:resolution_eta:1.5 -P SPH:CFL_condition:0.2"
../../../swift --hydro --cosmology --self-gravity --threads=8 $SWIFT_ARGS zeldovichPancake.yml 2>&1 | tee output.log

# Plot the result (one figure per snapshot)
for f in zeldovichPancake_[0-9][0-9][0-9][0-9].hdf5
do
    python3 plotSolution.py ${f:17:4}
done
