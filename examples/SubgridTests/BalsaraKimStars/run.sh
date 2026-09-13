#!/bin/bash

# Generate the initial conditions if not present
if [ ! -e glassCube_64.hdf5 ]
then
        echo "Fetching initial conditions for the Balsara-Kim test"
        wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/glassCube_64.hdf5
fi

if [ ! -e BalsaraKim.hdf5 ]
then
        echo "creating initial conditions"
        python3 makeIC.py
fi

# Get the cooling tables
if [ ! -e coolingtables ]
then
        echo "Fetchin cooling tables"
        wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/CoolingTables/EAGLE/coolingtables.tar.gz
        tar -xf coolingtables.tar.gz
fi

# Get the yield tables
if [ ! -e yieldtables ]
then
	echo "Fetching yield tables"
	wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/YieldTables/EAGLE/yieldtables.tar.gz
	tar -xf yieldtables.tar.gz
fi

# Get the photometry tables
if [ ! -e photometry ]
then
	echo "Fetching photometry tables"
	wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/YieldTables/EAGLE/photometry.tar.gz
	tar -xf photometry.tar.gz
fi

# Run SWIFT
../../../swift --hydro --external-gravity --cooling --stars --feedback --limiter --sync --threads=5 balsarakim.yml |& tee output.log
