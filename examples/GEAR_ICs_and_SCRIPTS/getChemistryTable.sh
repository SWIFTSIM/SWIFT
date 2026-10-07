#!/bin/bash

# Define target paths
DEST_DIR="./"

# Default state: Download with winds
WITH_WINDS=false

# Print usage instructions
usage() {
    echo "Usage: $0 [-n] [-h]"
    echo "  --with-winds     Download feedback tables with stellar winds"
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

    # TODO: the Data/SW groups of POPII.hdf5 and POPIII_PISNe.hdf5 must be
    # indexed by absolute Z, not Z/Zsun (see
    # stellar_evolution_compute_preSN_properties), and SWIFT stops at start-up
    # if a Data/SW group lacks metallicity_convention. Upload the regenerated
    # tables and update these URLs.
    wget -P "$DEST_DIR" https://obswww.unige.ch/~revazy/DATA/Swift/PreSNeTables/POPII.hdf5
    wget -P "$DEST_DIR" https://obswww.unige.ch/~revazy/DATA/Swift/PreSNeTables/POPIII_PISNe.hdf5

else
    echo "========================================"
    echo "Downloading feedback tables WITHOUT stellar winds..."
    echo "Source: Cosma Durham web storage"
    echo "========================================"

    wget -P "$DEST_DIR" https://virgodb.cosma.dur.ac.uk/swift-webstorage/FeedbackTables/POPIIsw.h5
fi

echo "Done."
