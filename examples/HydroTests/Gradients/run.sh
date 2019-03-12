#! /bin/bash

python makeICs.py stretched
../../swift --hydro --threads=2 gradientsStretched.yml
python plot.py gradients_stretched_0001.hdf5 stretched

python makeICs.py cartesian
../../swift --hydro --threads=2 gradientsCartesian.yml
python plot.py gradients_cartesian_0001.hdf5 cartesian

python makeICs.py random
../../swift --hydro --threads=2 gradientsRandom.yml
python plot.py gradients_random_0001.hdf5 random
