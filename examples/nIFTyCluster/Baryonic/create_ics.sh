#!/bin/bash

set -eo pipefail

# Convert the original (Gadget binary) nIFTy ICs to nifty.hdf5.
# Usage: ./create_ics.sh /path/to/IC_CLUSTER_00019
if [ $# -ne 1 ]
then
    echo "Usage: $0 /path/to/IC_CLUSTER_00019"
    exit 1
fi

bash setup.sh

source pygadenv/bin/activate

python3 convert_to_swift.py "$1"

deactivate

bash clean.sh
