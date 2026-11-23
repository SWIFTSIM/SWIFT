#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e cosmicRayBox.hdf5 ]
then
    echo "Generating the initial conditions of the cosmic-ray box..."
    python3 makeIC.py
fi

# Run SWIFT with the constant rate, then with the rate stored in each particle
../../../swift --hydro --cooling --threads=4 cosmicRayBox.yml
../../../swift --hydro --cooling --threads=4 cosmicRayBoxPerParticle.yml

# Compare with a one-zone run of pygeajay
python3 checkEquilibrium.py cosmicRayBox &&
    python3 checkEquilibrium.py cosmicRayBoxPerParticle --per-particle
