#!/bin/bash

# Generate the initial conditions if they are not present.
echo "Generating initial conditions for the cooling box example..."

python makeIC.py 10

../swift -s -t 1 coolingBox.yml -C 2>&1 | tee output.log

python energy_plot.py 0
