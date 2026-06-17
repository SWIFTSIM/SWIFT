#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e plummer.hdf5 ]
then
    echo "Generating initial conditions for the AGAMA potential example..."
    ic_plummer   -N 10000 --Mtot 1e-2 -a 0.1 --Rmax 20 -t swift -o plummer.hdf5
fi

../../../swift --external-gravity --threads=1 params.yml 2>&1 | tee output.log

python3 ./plotdensity.py snap/snapshot_0000.hdf5 snap/snapshot_0010.hdf5
python3 ./plotenergy.py  statistics.txt --x Time --y Rel_Tot_Energy -o energy.png

