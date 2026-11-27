#!/bin/bash

set -e
set -o pipefail

if [ ! -e zoom_uniform_dm_gravity_holes.hdf5 ]; then
	echo "Generating initial conditions for UniformDMGravityWithHoles..."
	python3 makeIC.py
fi

rm -f zoom_uniform_dm_gravity_holes_*.hdf5 statistics.txt output.log

../../../swift --cosmology --self-gravity --zoom --threads=8 -n 16 zoom_uniform_dm_gravity_holes.yml 2>&1 | tee output.log
