#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e dwarf_galaxy.hdf5 ]
then
    echo "Fetching initial conditions for the dwarf galaxy example..."
    ./getIC.sh
fi

../../swift --feedback --limiter --sync --self-gravity --hydro --stars --cooling --star-formation --threads=8 $@ dwarf_galaxy.yml 2>&1 | tee output.log

