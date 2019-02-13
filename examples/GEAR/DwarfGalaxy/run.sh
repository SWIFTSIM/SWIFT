#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e dwarf_galaxy.hdf5 ]
then
    echo "Fetching initial conditions for the dwarf galaxy example..."
    ./getIC.sh
fi

../../swift --feedback --self-gravity --hydro --stars --threads=8 $@ dwarf_galaxy.yml 2>&1 | tee output.log

