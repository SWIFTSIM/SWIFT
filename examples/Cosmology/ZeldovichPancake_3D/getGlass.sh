#!/bin/bash

set -eo pipefail

# Glass file for makeIC_glass.py (the lattice version makeIC.py needs none)
wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/glassCube_32.hdf5
