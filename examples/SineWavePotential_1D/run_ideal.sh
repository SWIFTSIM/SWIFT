#!/bin/bash

if [ ! -e sineWavePotential.hdf5 ]
then
  echo "Generating initial conditions for the 1D SineWavePotential example..."
  python makeIC_ideal.py
fi

../swift -g -s -t 2 sineWavePotential.yml 2>&1 | tee output.log

for f in sineWavePotential_*.hdf5
do
  python plotSolution_ideal.py $f
done
