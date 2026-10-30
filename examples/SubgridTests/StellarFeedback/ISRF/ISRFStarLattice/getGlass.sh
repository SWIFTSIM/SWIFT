#!/bin/bash
# Glass cube of n^3 particles, n=32.
n=${1:-32}
wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/glassCube_${n}.hdf5
