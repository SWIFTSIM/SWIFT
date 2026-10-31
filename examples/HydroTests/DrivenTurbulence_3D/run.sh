#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e drivenTurbulence.hdf5 ]
then
    echo "Generating initial conditions for the driven-turbulence example..."
    python3 makeIC.py 64
fi

# Run SWIFT (configure with --with-forcing=driven-turbulence --with-adiabatic-index=1.001)
../../../swift --hydro --limiter --threads=4 drivenTurbulence.yml 2>&1 | tee output.log

# Velocity power spectrum and Mach-number history against Bauer & Springel (2012)
python3 plotSolution.py
