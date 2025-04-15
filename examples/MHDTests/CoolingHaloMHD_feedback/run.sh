#!/bin/bash

# Generate the initial conditions if they are not present.
echo "Generating initial conditions for the MHD cloud"
python3 makeIC_MWsize.py 450000 

# Run SWIFT
../../../swift --threads=16 --feedback --external-gravity --self-gravity --stars --star-formation --cooling --hydro --limiter --sync cooling_halo.yml 2>&1 | tee output.log
