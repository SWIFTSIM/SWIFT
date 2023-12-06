#!/bin/bash

# make run.sh fail if a subcommand fails
set -e
set -o pipefail

if [ ! -e glassPlane_64.hdf5 ]
then
    echo "Fetching initial glass file ..."
    ./getGlass.sh
fi

if [ ! -f 'advect_metals.hdf5' ]; then
    echo "Generating ICs"
    python3 makeIC.py
fi

# Run SWIFT (must be compiled with --with-chemistry=EAGLE)
../../../swift \
    --hydro --threads=4 advect_metals.yml 2>&1 | tee output.log

python3 runSanityChecksAdvectedMetals.py
python3 plotAdvectedMetals.py
