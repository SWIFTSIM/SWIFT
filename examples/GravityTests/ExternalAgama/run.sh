#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e plummer.hdf5 ]
then
    echo "Generating initial conditions for the AGAMA potential example..."
    python3 plummerIC.py
fi

../../../swift --external-gravity --threads=1 params.yml 2>&1 | tee output.log

python3 ./plotdensity.py snap/snapshot_0000.hdf5 snap/snapshot_0010.hdf5
python3 ./plotenergy.py  statistics.txt --x Time --y Rel_Tot_Energy -o energy.png

