#!/bin/bash -l

# Ensure that the right modules are loaded, MAKE SURE to load python3 as well
module purge
# module load ...
# module load ...
# module load ...

# Create ICs
python3 write_ICs.py

../../../swift --hydro --temperature --threads=16 --limiter --sync --pin  idealized_jet.yml

