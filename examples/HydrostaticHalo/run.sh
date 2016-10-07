#!/bin/bash

# Generate the initial conditions if they are not present.
echo "Generating initial conditions for the isothermal potential box example..."
python makeIC.py 10000 

../swift -g -s -t 16 hydrostatic.yml 2>&1 | tee output.log

python radial_profile.py 10

python internal_energy_profile.py 10

python test_energy_conservation.py
