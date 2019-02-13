#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e initial_conditions.hdf5 ]
then
    echo "Generating initial conditions for the keplerian ring example..."
    echo "Please consider choosing your own options before continuing..."
    python3 makeIC.py
fi

rm -rf keplerian_ring_*.hdf5
../../swift --external-gravity --hydro --threads=1 --verbose=1 keplerian_ring.yml 2>&1 | tee output.log

echo
echo
echo "This test is expected to fail, please read the README.md file"
echo
echo
