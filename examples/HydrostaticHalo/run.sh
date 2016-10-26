#!/bin/bash

# Generate the initial conditions if they are not present.
echo "Generating initial conditions for the isothermal potential box example..."
python makeIC.py 10000 

../swift -g -s -t 16 hydrostatic.yml 2>&1 | tee output.log

python density_profile.py 2. 200 100

python internal_energy_profile.py 2. 200 100

python velocity_profile.py 2. 200 100

python test_energy_conservation.py
