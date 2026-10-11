#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e initial_conditions.hdf5 ]
then
    echo "Generating initial conditions for the Cullen & Dehnen (2010) Keplerian ring..."
    python3 makeIC.py
fi

# Run SWIFT (2D code, softened point-mass potential:
# ./configure --with-hydro-dimension=2 --with-ext-potential=point-mass-softened)
rm -rf keplerian_ring_*.hdf5
../../../swift --external-gravity --hydro --threads=4 keplerian_ring.yml 2>&1 | tee output.log
