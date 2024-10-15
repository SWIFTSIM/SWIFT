#!/bin/bash

# make run.sh fail if a subcommand fails
set -e
set -o pipefail

if [ ! -e glassPlane_64.hdf5 ]
then
    echo "Fetching initial glass file ..."
    ./getGlass.sh
fi

# Always generate ICs in case ELEMENT_COUNT has been changed
echo "Generating ICs"
python3 makeIC.py

# Run SWIFT (must be compiled with --with-chemistry=GEAR_X or --with-chemistry=AGORA)
../../../swift \
    --hydro --threads=4 advect_metals.yml 2>&1 | tee output.log

python3 runSanityChecksAdvectedMetals.py
python3 plotAdvectedMetals.py
