#!/bin/bash

 # Generate the initial conditions.
echo "Generating initial conditions for the Simple Orbits example..."
python makeIC.py

# Run SWIFT
../../swift --logger --external-gravity --threads=1 simple_orbits.yml 2>&1 | tee output.log

# Plot the solution
python3 plotSolution.py
