#!/bin/bash

# make run.sh fail if a subcommand fails
set -e
set -o pipefail

if [ ! -e glassCube_64.hdf5 ]
then
    echo "Fetching initial glass file for Strömgen Sphere 3D example ..."
    ./getGlass.sh
fi

if [ ! -f 'stromgrenSphere-3D.hdf5' ]; then
    echo "Generating ICs"
    python3 makeIC.py
fi

# Run SWIFT with RT
../../../swift \
    --hydro --threads=4 --stars --external-gravity \
    --feedback --radiation \
    stromgrenSphere-3D.yml 2>&1 | tee output.log

# option with mpi
# mpirun -np 2 ../../../swift_mpi --hydro --threads=14 --stars --external-gravity --feedback --radiation stromgrenSphere-3D.yml 2>&1 | tee output.log

# Plot the Stromgren 3D checks.
python3 ./plotSolution.py
