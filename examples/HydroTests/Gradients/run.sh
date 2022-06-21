#! /bin/bash

python3 makeICs.py stretched
../../swift --hydro --threads=2 gradientsStretched.yml
python3 plot.py gradients_stretched_0001.hdf5 stretched

python3 makeICs.py cartesian
../../swift --hydro --threads=2 gradientsCartesian.yml
python3 plot.py gradients_cartesian_0001.hdf5 cartesian

python3 makeICs.py random
../../swift --hydro --threads=2 gradientsRandom.yml
python3 plot.py gradients_random_0001.hdf5 random
