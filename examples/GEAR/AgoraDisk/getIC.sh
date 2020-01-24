#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "You need to provide the resolution (e.g. ./getIC.sh low)."
    echo "The possible options are low, med and high."
    exit
fi

wget http://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/AgoraDisk/$1.hdf5
