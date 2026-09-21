#!/bin/bash

# Default state (no flags given): download the standard Cloudy tables
# (HM2012 + high density). Passing any flag selects exactly the
# requested table(s) instead.
HM2012=false
HIGH_DENSITY=false
SHIELDED=false
ANY_FLAG=false

usage() {
    echo "Usage: $0 [--hm2012] [--high-density] [--shielded] [-h]"
    echo "  --hm2012        Download the standard Cloudy table"
    echo "  --high-density  Download the high density Cloudy table"
    echo "  --shielded      Download the GEAR shielded Cloudy table"
    echo "  -h, --help      Display this help menu"
    echo "  No flags: downloads --hm2012 and --high-density."
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
	--hm2012)
	    HM2012=true
	    ANY_FLAG=true
	    shift
	    ;;
	--high-density)
	    HIGH_DENSITY=true
	    ANY_FLAG=true
	    shift
	    ;;
	--shielded)
	    SHIELDED=true
	    ANY_FLAG=true
	    shift
	    ;;
	-h|--help)
	    usage
	    ;;
	*)
	    echo "Unknown option: $1"
	    usage
	    ;;
    esac
done

if [ "$ANY_FLAG" = false ]; then
    HM2012=true
    HIGH_DENSITY=true
fi

if [ "$HM2012" = true ]; then
    wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/CoolingTables/CloudyData_UVB=HM2012.h5
fi

if [ "$HIGH_DENSITY" = true ]; then
    wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/CoolingTables/CloudyData_UVB=HM2012_high_density.h5
fi

if [ "$SHIELDED" = true ]; then
    wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/CoolingTables/GEAR/CloudyData_UVB=HM2012_shielded.h5
fi
