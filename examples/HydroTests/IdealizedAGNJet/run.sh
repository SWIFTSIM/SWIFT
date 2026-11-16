#!/bin/bash -l

# Ensure that the right modules are loaded, MAKE SURE to load python3 as well
#module purge

set -eo pipefail

# Create ICs if they are not present (remove ICs.hdf5 after changing write_ICs.py)
if [ ! -e ICs.hdf5 ]
then
    python3 write_ICs.py
fi

../../../swift --hydro --temperature --threads=16 --limiter --sync --pin  idealized_jet.yml

cd analysis_scripts

python3 energetics.py
python3 lobe_dimensions.py
python3 phase_space_diagram.py
python3 temperature_slices.py
