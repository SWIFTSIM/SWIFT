
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

# Run SWIFT
../swift -s -C -t 1 coolingBox.yml

# Check energy conservation and cooling rate
python energy_plot.py
