#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e square_equal_spacing.hdf5 ]
then
    echo "Generating initial conditions for the 3D square test (equal spacing)..."
    python3 makeICEqualSpacing.py
fi
if [ ! -e square_equal_mass.hdf5 ]
then
    echo "Generating initial conditions for the 3D square test (equal mass)..."
    python3 makeICEqualMass.py
fi

# Run SWIFT
../../../swift --hydro --threads=4 square_equal_spacing.yml 2>&1 | tee output_equal_spacing.log
../../../swift --hydro --threads=4 square_equal_mass.yml 2>&1 | tee output_equal_mass.log

# Plot the solutions
python3 ./plotSolution.py "equal_spacing" 20
python3 ./plotSolution.py "equal_mass" 20
