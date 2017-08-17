#! /bin/bash

python makeICs.py stretched
../swift -s -t 2 gradientsStretched.yml
python plot.py gradients_stretched_0001.hdf5 stretched

python makeICs.py cartesian
../swift -s -t 2 gradientsCartesian.yml
python plot.py gradients_cartesian_0001.hdf5 cartesian

python makeICs.py random
../swift -s -t 2 gradientsRandom.yml
python plot.py gradients_random_0001.hdf5 random
