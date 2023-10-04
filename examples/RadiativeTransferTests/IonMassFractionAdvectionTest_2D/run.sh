#!/bin/bash

# make run.sh fail if a subcommand fails
set -e
set -o pipefail

if [ ! -e glassPlane_64.hdf5 ]
then
    echo "Fetching initial glass file ..."
    ./getGlass.sh
fi

if [ ! -f 'advect_ions.hdf5' ]; then
    echo "Generating ICs"
    python3 makeIC.py
fi

# Run SWIFT with RT
../../../swift \
    --hydro --threads=4 --stars --external-gravity \
    --feedback --radiation \
    advect_ions.yml 2>&1 | tee output.log

python3 plotIonization.py
python3 testIonization.py
