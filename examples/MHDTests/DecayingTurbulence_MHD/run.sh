#!/bin/bash

# Generate the initial conditions file if not present.
if [ ! -f ./DecayingTurbulence.hdf5 ]
then
    echo "Generating a BCC unit cell, copies of which are to be stacked to generate the left state of the Brio & Wu shock tube ..."
    python3 makeIC.py
fi

# Run the example with SWIFT 
../../../swift --hydro --threads=4 DecayingTurbulence.yml 2>&1 | tee output.log

# Plot the calculated solution at time t=0.1
python3 plotSpectrum.py DecayingTurbulence_0011.hdf5 DecayingTurbulence_0011.png
