#!/bin/bash

if [ ! -e particleSplitting.hdf5 ]
then
    echo "Making initial conditions for the ParticleSplitting example"
    python3 makeIC.py
fi


# A very simple invocation for this one...
../../../swift --hydro -t 4 particle_splitting.yml

python3 plotSolution.py
