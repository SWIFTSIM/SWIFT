#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e FastRotor_LR.hdf5 ]
then
    echo "Fetching glass files for generating ICs..."
    ./getGlass.sh
    echo "Generating the Bfield"
    python ./makeIC_VP.py
fi

cd VeP
./run.sh
cd ../ODI
./run.sh
cd ../FDI
./run.sh


