#!/bin/bash

# Generate the initial conditions if they are not present.
echo "Generating initial conditions for the uniform box example..."

python makeIC.py 10

../swift -s -C -t 16 uniformBox.yml

python energy_plot.py 0
