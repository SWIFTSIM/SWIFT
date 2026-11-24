#!/bin/bash
# Glass cube of n^3 particles, n=64.
n=${1:-64}
wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/glassCube_${n}.hdf5
