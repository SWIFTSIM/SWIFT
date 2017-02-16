#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e Isothermal.hdf5 ]
then
    echo "Generating initial conditions for the isothermal potential box example..."
    python makeIC.py 1000 1
fi

rm -rf Isothermal_*.hdf5
../swift -g -t 1 isothermal.yml 2>&1 | tee output.log

python energy_plot.py
