#!/bin/bash
#
# Put a GEAR yields table carrying the radiation data into the current
# directory.
#
# Usage: getRadiationTable.sh [table.hdf5]
#
# Photoionization, radiation pressure and the interstellar radiation field
# read the stellar photon rates and luminosities from a "Data/Radiation"
# group. No public host serves such a table yet: the files
# getChemistryTable.sh downloads predate it and carry no radiation data at
# all, so that script cannot supply one. The table is generated locally
# with pychem instead, and this script only resolves it into the example
# directory.
#
# Resolution order:
#   1. $GEAR_RADIATION_TABLE, if set, is the source file.
#   2. $HOME/programs/pychem/<name>, the default pychem output location.
# Otherwise the script explains how to generate one and fails.

set -eu

table="${1:-PopII_parsec_spectral.hdf5}"

if [ -e "$table" ]; then
    exit 0
fi

source_file="${GEAR_RADIATION_TABLE:-$HOME/programs/pychem/$table}"

if [ -e "$source_file" ]; then
    echo "Using the radiation yields table at '$source_file'."
    cp "$source_file" "$table"
    exit 0
fi

cat >&2 <<EOF

ERROR: no radiation yields table found for this example.

Looked for:
  ./$table
  $source_file

Generate one with pychem's pychem_generate_hdf5_parameters, using a
spectral parameter file (the table's Data/Radiation group then carries
qh_source='spectral'), then either copy it here under the name above or
point GEAR_RADIATION_TABLE at it:

  GEAR_RADIATION_TABLE=/path/to/table.hdf5 ./run.sh

The tables served by the public hosts (see getChemistryTable.sh) carry no
Data/Radiation group and cannot drive this example.

EOF
exit 1
