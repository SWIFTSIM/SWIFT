#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e EAGLE_ICs_12.hdf5 ]
then
    echo "Fetching initial conditions for the EAGLE 12Mpc example..."
    ./getIC.sh
fi

mpirun -np 4 ../swift_mpi -s -c -C -t 16 -n 100 eagle_12.yml 2>&1 | tee output.log

