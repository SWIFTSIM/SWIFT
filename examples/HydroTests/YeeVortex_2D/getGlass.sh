#!/bin/bash

set -eo pipefail
# 2D glass files: 128^2 for the default run, 32^2 and 64^2 for the lower resolutions of the convergence study.
for n in 32 64 128; do
    [ -e glassPlane_$n.hdf5 ] || wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/glassPlane_$n.hdf5
done
