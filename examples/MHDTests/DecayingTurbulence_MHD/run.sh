#!/bin/bash

# Generate the initial conditions file if not present.
if [ ! -f ./DecayingTurbulence.hdf5 ]
then
    echo "Generating a unit with two excited modes: one large scale velocity and small scale MF"
    python3 makeIC.py
fi

# Run the example with SWIFT 
../../../swift --hydro --threads=4 DecayingTurbulence.yml 2>&1 | tee output.log

# Plot the calculated solution at time t=0.1
python3 plotSpectrum.py DecayingTurbulence_0011.hdf5 DecayingTurbulence_0020.png
