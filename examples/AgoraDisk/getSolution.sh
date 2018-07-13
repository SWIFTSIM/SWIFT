#!/bin/bash

OPTIND=1

with_cooling=0

function show_help {
    echo "Valid options are:"
    echo "\t -h \t Show this help"
    echo "\t -C \t Download solution with cooling"
}

while getopts "h?C" opt; do
    case "$opt" in
	h|\?)
	    show_help
	    exit
	    ;;
	C)
	    with_cooling=1
	    ;;
    esac
done

# cleanup work space
rm snapshot_0000.hdf5 snapshot_0500.hdf5

if [ $with_cooling -eq 1 ]; then
    wget https://obswww.unige.ch/~lhausamm/swift/IC/AgoraDisk/Gear/with_cooling/snapshot_0000.hdf5
    wget https://obswww.unige.ch/~lhausamm/swift/IC/AgoraDisk/Gear/with_cooling/snapshot_0500.hdf5
else
    wget https://obswww.unige.ch/~lhausamm/swift/IC/AgoraDisk/Gear/without_cooling/snapshot_0000.hdf5
    wget https://obswww.unige.ch/~lhausamm/swift/IC/AgoraDisk/Gear/without_cooling/snapshot_0500.hdf5
fi


