#!/bin/bash
# Glass cube size must match the requested resolution level: N = n^3
# particles, n = 2**level (see makeIC.py). Defaults to 32 (level=5).
n=${1:-32}
wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/glassCube_${n}.hdf5
