#!/bin/bash

# Generate the initial conditions if they are not present.
echo "Generating initial conditions for the MHD isothermal potential box example..."
python3 makeIC_dwarfsize.py 1000000 

# Run SWIFT with external potential, SPH and cooling
../../../swift --external-gravity --hydro --cooling --threads=70 cooling_halo_dwarfsize.yml 2>&1 | tee output.log
