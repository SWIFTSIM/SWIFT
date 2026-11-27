#!/bin/bash

set -e
set -o pipefail

if [ ! -e boundary_offset_uni_dm_grav.hdf5 ]; then
	echo "Generating initial conditions for BoundaryOffsetUniDMGrav..."
	python3 makeIC.py
fi

rm -f boundary_offset_uni_dm_grav_*.hdf5 statistics.txt output.log

../../../swift --cosmology --self-gravity --zoom --threads=8 -n 16 boundary_offset_uni_dm_grav.yml 2>&1 | tee output.log
