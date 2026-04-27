#!/bin/bash

set -e
set -o pipefail

if [ ! -e zoom_uniform_dm_hydro.hdf5 ]; then
	echo "Generating initial conditions for UniformDMHydro..."
	python3 makeIC.py
fi

rm -f zoom_uniform_dm_hydro_*.hdf5 statistics.txt output.log

../../../swift --cosmology --hydro --self-gravity --zoom --threads=8 -n 16 zoom_uniform_dm_hydro.yml 2>&1 | tee output.log
