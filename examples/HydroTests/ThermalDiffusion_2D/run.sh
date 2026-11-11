#!/bin/bash

set -eo pipefail

# Generate the initial conditions if they are not present.
if [ ! -e glassPlane_128.hdf5 ]
then
    echo "Fetching initial glass file for the thermal diffusion example..."
    ./getGlass.sh
fi
if [ ! -e thermalDiffusion.hdf5 ]
then
    echo "Generating initial conditions for the thermal diffusion example..."
    python3 makeIC.py
fi

# Run SWIFT (configure with --with-hydro-dimension=2)
../../../swift --hydro --threads=4 thermalDiffusion.yml 2>&1 | tee output.log

# Profiles and the effective conductivity
python3 plotSolution.py
