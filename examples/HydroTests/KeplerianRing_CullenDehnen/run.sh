#!/bin/bash

set -eo pipefail

# Generate the initial conditions if they are not present.
if [ ! -e initial_conditions.hdf5 ]
then
    echo "Generating initial conditions for the Cullen & Dehnen (2010) Keplerian ring..."
    python3 makeIC.py
fi

# Run SWIFT (2D code, softened point-mass potential:
# ./configure --with-hydro-dimension=2 --with-ext-potential=point-mass-softened)
# SWIFT_ARGS passes extra arguments, e.g. MAGMA's recommended
# SWIFT_ARGS="-P SPH:resolution_eta:1.5 -P SPH:CFL_condition:0.2"
../../../swift --external-gravity --hydro --threads=4 $SWIFT_ARGS keplerian_ring.yml 2>&1 | tee output.log
