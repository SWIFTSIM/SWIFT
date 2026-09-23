#!/bin/bash

# Define target paths
DEST_DIR="./"

# Default state: Download with winds
WITH_WINDS=false

# Default state: do not fetch the radiation tables
WITH_RADIATION=false

# Print usage instructions
usage() {
    echo "Usage: $0 [-n] [-h]"
    echo "  --with-winds     Download feedback tables with stellar winds"
    echo "  --with-radiation Download the spectral tables carrying Data/Radiation"
    echo "  -h, --help      Display this help menu"
    exit 1
}

# Parse command line flags
while [[ $# -gt 0 ]]; do
    case "$1" in
	--with-winds)
	    WITH_WINDS=true
	    shift
	    ;;
	--with-radiation)
	    WITH_RADIATION=true
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

# Execute download based on configuration flag
if [ "$WITH_WINDS" = true ]; then
    echo "========================================"
    echo "Downloading feedback tables WITH stellar winds..."
    echo "Source: UniGe Astro servers"
    echo "========================================"

    wget -P "$DEST_DIR" https://obswww.unige.ch/~revazy/DATA/Swift/PreSNeTables/POPII.hdf5
    wget -P "$DEST_DIR" https://obswww.unige.ch/~revazy/DATA/Swift/PreSNeTables/POPIII_PISNe.hdf5

else
    echo "========================================"
    echo "Downloading feedback tables WITHOUT stellar winds..."
    echo "Source: Cosma Durham web storage"
    echo "========================================"

    wget -P "$DEST_DIR" https://virgodb.cosma.dur.ac.uk/swift-webstorage/FeedbackTables/POPIIsw.h5
fi

if [ "$WITH_RADIATION" = true ]; then
    echo "========================================"
    echo "Downloading the spectral tables WITH radiation data..."
    echo "========================================"

    # One resolver for these two, so the share ids and their checksums
    # live in a single place.
    "$(dirname "$0")"/getRadiationTable.sh PopII_parsec_spectral.hdf5
    "$(dirname "$0")"/getRadiationTable.sh PopIII_parsec_spectral.hdf5
fi

echo "Done."
echo
echo "Note: GEAR photoionization, radiation pressure and the interstellar"
echo "radiation field need a table carrying a 'Data/Radiation' group. The"
echo "public-host tables above do not. Pass --with-radiation for tables"
echo "that do, or generate one with pychem's"
echo "pychem_generate_hdf5_parameters. Check any table with"
echo "  ./checkRadiationTable.sh <table.h5> [--with-isrf]"
