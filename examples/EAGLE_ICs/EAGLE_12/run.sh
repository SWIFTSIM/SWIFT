#!/bin/bash

set -eo pipefail

# Generate the initial conditions if they are not present.
if [ ! -e EAGLE_L0012N0188_ICs.hdf5 ]
then
    echo "Fetching initial conditions for the EAGLE 12Mpc example..."
    ./getIC.sh
fi

# Grab the cooling, yield, and photometry tables if they are not present.
if [ ! -e yieldtables ]
then
    echo "Fetching EAGLE yield tables..."
    ../getEagleYieldTable.sh
fi

if [ ! -e UV_dust1_CR1_G1_shield1.hdf5 ]
then
    echo "Fetching EAGLE-XL cooling tables..."
    ../getPS2020CoolingTables.sh
fi

if [ ! -e photometry ]
then
    echo "Fetching EAGLE photometry tables..."
    ../getEaglePhotometryTable.sh
fi

# The run-time options are broken down by line as:
# Cosmology and all the EAGLE model options (--eagle is equivalent to
#   --hydro --limiter --sync --self-gravity --stars --star-formation
#   --cooling --feedback --black-holes --fof; fof is needed for the BH seeding)
# Threading options - run with threads and pinning (latter not required but improves performance)
# The corresponding parameter file for this run

../../../swift \
    --cosmology --eagle \
    --threads=16 --pin \
    eagle_12.yml

