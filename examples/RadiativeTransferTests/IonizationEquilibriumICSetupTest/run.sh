#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

if [ ! -f ./ionization_equilibrium_test.hdf5 ]; then
    echo "creating ICs"
    python3 makeIC.py
fi

# Run SWIFT with RT
../../swift \
    --hydro \
    --threads=1 \
    --verbose=0  \
    --radiation \
    --stars \
    --feedback \
    --external-gravity \
    --steps 1 \
    ./ion_equil.yml 2>&1 | tee output.log


if [ ! -f "IonizationEquilibriumICSetupTestReference.txt" ]; then
    ./getReference.sh
fi
python3 plotSolution.py
