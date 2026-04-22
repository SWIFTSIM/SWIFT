#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e EAGLE_ICs_6.hdf5 ]
then
    echo "Fetching initial conditions for the EAGLE 6Mpc example..."
    ./getIC.sh
fi

if [ ! -e EAGLE_ICs_6_SIDM.hdf5 ]
then
    python3 redo_ICs.py EAGLE_ICs_6.hdf5 EAGLE_ICs_6_SIDM.hdf5 
fi

echo "Running with SIDM..."
../../../../swift -v 1 --sidm --self-gravity --threads=4 -n 1 -y 1 eagle_6_sidm.yml | tee output_sidm.log

echo "Running with gas..."
../../../../swift -v 1 --hydro --threads=4 -n 1 -y 1 eagle_6.yml | tee output.log

echo "Comparing gas and SIDM densities"
python3 compare_gas_SIDM_densities.py --sidm eagle_sidm_0000.hdf5 --gas eagle_0000.hdf5