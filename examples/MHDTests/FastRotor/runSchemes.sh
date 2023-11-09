#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e FastRotor_LR.hdf5 ]
then
    echo "Fetching Glass Files..."
    ./getGlass.sh
    echo "Generating the ICs"
    python ../makeIC_VP.py 
fi

cd VeP
./run.sh &
cd ../ODI
./run.sh &
cd ../FDI
./run.sh
cd ..
echo "RUNING .... "
