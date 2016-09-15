#!/bin/bash

# Generate the initial conditions if they are not present.
echo "Generating initial conditions for the isothermal potential box example..."
python makeIC.py 1000 

../swift -g -s -t 16 hydrostatic.yml 2>&1 | tee output.log
