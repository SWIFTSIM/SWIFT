#!/bin/bash

set -e
set -o pipefail

if [ ! -e offset_uni_dm_gravity.hdf5 ]; then
    echo "Generating initial conditions for OffsetUniDMGravity..."
    python3 makeIC.py
fi

rm -f offset_uni_dm_gravity_*.hdf5 statistics.txt output.log

../../../swift --cosmology --self-gravity --zoom --threads=8 -n 16 -v 1 offset_uni_dm_gravity.yml 2>&1 | tee output.log
