#!/bin/bash

set -e
set -o pipefail

if [ ! -e zoom_uniform_dm_eagle.hdf5 ]; then
	echo "Generating initial conditions for UniformDMEAGLE..."
	python3 makeIC.py
fi

if [ ! -e yieldtables ]; then
	echo "Fetching EAGLE yield tables..."
	../../EAGLE_ICs/getEagleYieldTable.sh
fi

if [ ! -e coolingtables ]; then
	echo "Fetching EAGLE cooling tables..."
	../../EAGLE_ICs/getEagleCoolingTable.sh
fi

if [ ! -e UV_dust1_CR1_G1_shield1.hdf5 ]; then
	echo "Fetching PS2020 cooling tables..."
	../../EAGLE_ICs/getPS2020CoolingTables.sh
fi

if [ ! -e photometry ]; then
	echo "Fetching EAGLE photometry tables..."
	../../EAGLE_ICs/getEaglePhotometryTable.sh
fi

rm -f zoom_uniform_dm_eagle_*.hdf5 statistics.txt output.log

../../../swift --cosmology --eagle --zoom --threads=4 -n 16 zoom_uniform_dm_eagle.yml 2>&1 | tee output.log
