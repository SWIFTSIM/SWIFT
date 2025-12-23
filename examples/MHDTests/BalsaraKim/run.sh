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

# Run SWIFT
../../../swift --hydro --cooling --threads=8 balsarakim.yml |& tee output.log
