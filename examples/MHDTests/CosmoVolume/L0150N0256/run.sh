#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e SWIFT_ICs_L0150N0256.hdf5 ]
then
    echo "Fetching initial conditions for the small cosmological volume example..."
    ./getIC.sh
fi

if [ ! -e SWIFT_MHD_ICs_L0150N0256.hdf5 ]
then
echo "Generating the Bfield"
    python3 ../makeIC.py SWIFT_ICs_L0150N0256.hdf5 SWIFT_MHD_ICs_L0150N0256.hdf5
fi

# Run SWIFT
../../../../swift --cosmology --hydro --self-gravity --fof --power --threads=16 swift_param_L0150N0256.yml 2>&1 | tee output.log

# Plot the temperature evolution
# python3 plotTempEvolution.py