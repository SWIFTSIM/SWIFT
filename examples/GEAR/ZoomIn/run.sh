#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e zoom_in.hdf5 ]
then
    echo "Fetching initial conditions for the zoom in example..."
    ./getIC.sh
fi

../../swift --cooling --feedback --cosmology  --limiter --sync --self-gravity --hydro --stars --star-formation --threads=8 zoom_in.yml 2>&1 | tee output.log

