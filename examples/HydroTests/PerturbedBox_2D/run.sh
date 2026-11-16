#!/bin/bash

# This example requires SWIFT to be configured with:
#   --with-hydro-dimension=2

set -eo pipefail

# Generate the initial conditions if they are not present.
if [ ! -e perturbedPlane.hdf5 ]
then
    echo "Generating initial conditions for the perturbed plane example..."
    python3 makeIC.py 100
fi

../../../swift --hydro --threads=4 perturbedPlane.yml 2>&1 | tee output.log
