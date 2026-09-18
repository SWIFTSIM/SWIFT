#!/bin/bash
set -e

if [ ! -e SWIFT_IC_TurbBox_N32.hdf5 ]
then
    echo "Fetching the ICs..."
    ./getIC.sh
fi

../../../swift --hydro --self-gravity --threads=8 params.yml 2>&1 | tee output.log

python3 metal_projection.py --log --vmin -5 --vmax -2.5 --metal total snap/snapshot_*.hdf5


