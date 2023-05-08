#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "You need to provide the resolution (e.g. ./getIC.sh LOW)."
    echo "The possible options are LOW and MED."
    exit
fi

wget http://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/AgoraDisk/snap.$1.hdf5
mv snap.$1.hdf5 agora_disk.hdf5

