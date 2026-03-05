#!/bin/bash

# Generate the initial conditions if they are not present.
echo "Generating initial conditions for the basic SIDM example..."
python3 makeICs.py

# Run SWIFT
../../../swift --sidm --self-gravity --verbose 1 params.yml 2>&1 | tee output.log
