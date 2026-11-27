#!/bin/bash

# This example requires SWIFT to be configured with:
#   --with-ext-potential=sine-wave --with-equation-of-state=isothermal-gas

set -eo pipefail

if [ ! -e sineWavePotential.hdf5 ]
then
  echo "Generating initial conditions for the 3D SineWavePotential example..."
  python3 makeIC.py
fi

../../../swift --external-gravity --hydro --threads=2 sineWavePotential.yml 2>&1 | tee output.log

for f in sineWavePotential_*.hdf5
do
  python3 plotSolution.py $f
done
