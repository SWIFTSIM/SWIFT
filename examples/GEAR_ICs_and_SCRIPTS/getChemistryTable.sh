#!/bin/bash

set -eu

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

# Fetch $2 (URL) into $DEST_DIR/$1, skipping a file already there (a local
# copy, e.g. the operator's own radiation-carrying table, is left
# untouched) and failing loudly on an empty or truncated download instead
# of leaving a file that only fails later at H5Fopen with no clue why.
fetch() {
    name="$1"
    url="$2"
    dest="$DEST_DIR$name"

    if [ -e "$dest" ]; then
	return 0
    fi

    tmp="$dest.part"
    rm -f "$tmp"
    if ! wget -O "$tmp" "$url"; then
	rm -f "$tmp"
	echo "$0: could not download '$name' from $url." >&2
	exit 1
    fi
    if [ ! -s "$tmp" ]; then
	rm -f "$tmp"
	echo "$0: '$name' downloaded from $url is empty." >&2
	exit 1
    fi
    mv "$tmp" "$dest"
}

# Execute download based on configuration flag
if [ "$WITH_WINDS" = true ]; then
    echo "========================================"
    echo "Downloading feedback tables WITH stellar winds..."
    echo "Source: UniGe Astro servers"
    echo "========================================"

    fetch POPII.hdf5 https://obswww.unige.ch/~revazy/DATA/Swift/PreSNeTables/POPII.hdf5
    fetch POPIII_PISNe.hdf5 https://obswww.unige.ch/~revazy/DATA/Swift/PreSNeTables/POPIII_PISNe.hdf5

else
    echo "========================================"
    echo "Downloading feedback tables WITHOUT stellar winds..."
    echo "Source: Cosma Durham web storage"
    echo "========================================"

    fetch POPIIsw.h5 https://virgodb.cosma.dur.ac.uk/swift-webstorage/FeedbackTables/POPIIsw.h5
fi

if [ "$WITH_RADIATION" = true ]; then
    echo "========================================"
    echo "Downloading the spectral tables WITH radiation data..."
    echo "========================================"

    # One resolver for these two, so the share ids and their checksums
    # live in a single place.
    "$(dirname "$0")"/getRadiationTable.sh PopII_parsec_spectral.hdf5 || exit 1
    "$(dirname "$0")"/getRadiationTable.sh PopIII_parsec_spectral.hdf5 || exit 1
fi

echo "Done."
echo
echo "Note: GEAR photoionization, radiation pressure and the interstellar"
echo "radiation field need a table carrying a 'Data/Radiation' group. The"
echo "public-host tables above do not. Pass --with-radiation for tables"
echo "that do, or generate one with pychem's"
echo "pychem_generate_hdf5_parameters. Check any table with"
echo "  ./checkRadiationTable.sh <table.h5> [--with-isrf]"
