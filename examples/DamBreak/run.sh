#!/bin/bash

set -eo pipefail

# Generate the initial conditions if they are not present.
if [ ! -e dam_break.hdf5 ]
then
    echo "Generating initial conditions for the 3D dam break test..."
    python3 makeIC.py --n_part 1000000
fi

# Run SWIFT
../../swift --hydro --external-gravity --threads=16 dam_break.yml 2>&1 | tee output_dam_break.log

# Plot the solutions
python3 ./plotSnapshots.py
