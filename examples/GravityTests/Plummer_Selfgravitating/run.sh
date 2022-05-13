#!/bin/bash

if [ ! -e plummer.hdf5 ]
then
    echo "Generating initial conditions for Plummer example..."
    python3 plummerIC.py -a 0.1 -N 65536
fi

# create output directory
if [ ! -e snap ]
then
  mkdir snap
fi

../../swift --self-gravity --threads=8 params.yml 2>&1 | tee output.log

echo "Plotting results..."
# If params.yml is left at default values, should produce 10 snapshots
python3 plotdensity.py snap/snapshot_*.hdf5
