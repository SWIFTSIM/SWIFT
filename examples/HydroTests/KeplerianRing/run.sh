#!/bin/bash

set -eo pipefail

# Generate the initial conditions if they are not present.
if [ ! -e initial_conditions.hdf5 ]
then
    echo "Generating initial conditions for the keplerian ring example..."
    echo "Please consider choosing your own options before continuing..."
    python3 makeIC.py
fi

# Run SWIFT (configure with --with-ext-potential=point-mass-softened)
../../../swift --external-gravity --hydro --threads=4 keplerian_ring.yml 2>&1 | tee output.log

# Surface density and particle positions
python3 plotSolution.py
