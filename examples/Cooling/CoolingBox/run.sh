
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

# Get the Grackle cooling table
if [ ! -e CloudyData_UVB=HM2012.h5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    ../getCoolingTable.sh
fi

# Run SWIFT
../../swift --hydro --cooling --threads=4 -n 1000 coolingBox.yml

# Check energy conservation and cooling rate
python plotEnergy.py
