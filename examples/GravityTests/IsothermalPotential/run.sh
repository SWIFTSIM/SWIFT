#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e Isothermal.hdf5 ]
then
    echo "Generating initial conditions for the isothermal potential box example..."
    python makeIC.py 1000 0
fi

rm -rf Isothermal_*.hdf5
../../swift --external-gravity --threads=1 isothermal.yml 2>&1 | tee output.log

python energy_plot.py
