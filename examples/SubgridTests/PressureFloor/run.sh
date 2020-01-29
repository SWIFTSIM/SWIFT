
#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e glassCube_32.hdf5 ]
then
    echo "Fetching initial glass file for the cooling box example..."
    ./getGlass.sh
fi
if [ ! -e coolingBox.hdf5 ]
then
    echo "Generating initial conditions for the cooling box example..."
    python makeIC.py
fi

rm pressureFloor_*
# Run SWIFT
../../swift --self-gravity --hydro --cooling --threads=8 pressureFloor.yml

# Check if the simulation collapsed
python plotDensity.py 80
