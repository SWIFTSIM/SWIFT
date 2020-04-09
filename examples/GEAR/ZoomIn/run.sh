#!/bin/bash

echo "Fetching initial conditions for the zoom in example..."
./getIC.sh

../../swift --cooling --feedback --cosmology  --limiter --sync --self-gravity --hydro --stars --star-formation --threads=8 zoom_in.yml 2>&1 | tee output.log

echo "Fetching the solution from GEAR..."
./getSolution.sh

# Convert SWIFT snapshot to GEAR
python3 cleanupSwift.py snapshot_1010.hdf5 swift_final.hdf5


# Make the plots
./make_image.py
