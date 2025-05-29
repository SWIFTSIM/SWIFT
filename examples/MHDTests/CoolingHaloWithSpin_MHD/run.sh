#!/bin/bash

# Generate the initial conditions if they are not present.
echo "Generating initial conditions for the MHD isothermal potential box example..."
python3 makeIC.py 93800 

# Run SWIFT with external potential, cooling
#../../../swift --external-gravity --self-gravity --hydro --cooling --stars --star-formation --threads=4 cooling_halo.yml 2>&1 | tee output.log

# Run SWIFT with external potential, cooling, self-gravity and star formation
#../../../swift --external-gravity --self-gravity --hydro --cooling --stars --star-formation --threads=4 cooling_halo.yml 2>&1 | tee output.log

# Run SWIFT with external potential, cooling, self-gravity, star formation and feedback
../../../swift --external-gravity --self-gravity --hydro --cooling --stars --star-formation --feedback --sync --threads=4 cooling_halo.yml 2>&1 | tee output.log
