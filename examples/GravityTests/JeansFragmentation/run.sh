#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e jeansbox.hdf5 ]
then
    echo "Generating initial conditions for the adiabatic Jeans fragmentation example..."
    python3 makeIC.py --lJ 0.5  -o jeansbox.hdf5
fi

# create output directory
if [ ! -e snap ]
then
  mkdir snap
fi

# Run for some sound crossing time
../../swift --hydro --self-gravity --threads=1 jeansfragmentation.yml 2>&1 | tee output.log
