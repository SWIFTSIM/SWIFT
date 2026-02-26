#!/bin/bash

set -e
set -o pipefail

if [ ! -e zoom_uniform_dm_gravity.hdf5 ]; then
	echo "Generating initial conditions for UniformDMGravity..."
	python3 makeIC.py
fi

rm -f zoom_uniform_dm_gravity_*.hdf5 statistics.txt output.log

../../../swift --cosmology --self-gravity --zoom --threads=4 -n 16 zoom_uniform_dm_gravity.yml 2>&1 | tee output.log
