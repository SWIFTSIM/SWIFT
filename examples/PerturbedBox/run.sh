#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e perturbedBox.hdf5 ]
then
    echo "Generating initial conditions for the perturbed box example..."
    python makeIC.py 50
fi

../swift -s -t 16 uniformBox.yml
