#!/bin/bash
#
# Put a GEAR yields table carrying the radiation data into the current
# directory.
#
# Usage: getRadiationTable.sh [table.hdf5]
#
# Photoionization, radiation pressure and the interstellar radiation field
# read the stellar photon rates and luminosities from a "Data/Radiation"
# group. The tables getChemistryTable.sh downloads from the public hosts
# predate that group and carry no radiation data at all.
#
# Resolution order:
#   1. the file already in this directory, if present;
#   2. $GEAR_RADIATION_TABLE, if set;
#   3. $HOME/programs/pychem/<name>, the default pychem output location;
#   4. the SWITCHdrive share below, whose content is verified against the
#      SHA-256 recorded here.
# Otherwise the script explains how to generate one and fails.

set -eu

table="${1:-PopII_parsec_spectral.hdf5}"

# Published spectral tables: pychem, PARSEC stellar evolution, spectral
# Q_H and L_PE/L_LW. The checksum is what makes a truncated download or a
# silently replaced share fail loudly instead of running.
case "$table" in
    PopII_parsec_spectral.hdf5)
	share_id="ydicQptiff7WspX"
	sha256="72f4447ad454525ce7035d312a9bb400a049d3b8b455ec9dadc9c1e84742e5fc"
	;;
    PopIII_parsec_spectral.hdf5)
	share_id="9D5yQEAa3NNg4F8"
	sha256="de2c8f28cb8f596796ca93cad2694e5c04eb6a19bb638ffed418a7f3859404ff"
	;;
    # Mass-only ("M"), blackbody-fits Q_H/L table. The HIIRegions and
    # RadiationPressure examples' star masses are calibrated against this
    # table's own Q_H (e.g. Starbench's 26.75 Msun reproduces Bisbas et
    # al. 2015's 1e49 photons/s), so those examples require this exact
    # file, not the spectral one above.
    radiation_fits_popII.hdf5)
	share_id="2bdJyajQjwHnEeB"
	sha256="3bddd6d06feefdd6efa4d6b65b61631b96ccf0817f145dfa8523a5c63f5772e5"
	;;
    *)
	share_id=""
	sha256=""
	;;
esac

if [ -e "$table" ]; then
    exit 0
fi

source_file="${GEAR_RADIATION_TABLE:-$HOME/programs/pychem/$table}"

if [ -e "$source_file" ]; then
    echo "Using the radiation yields table at '$source_file'."
    cp "$source_file" "$table"
    exit 0
fi

if [ -n "$share_id" ]; then
    url="https://drive.switch.ch/index.php/s/$share_id/download"
    echo "Downloading the radiation yields table '$table'."
    tmp="$table.part"
    rm -f "$tmp"
    if command -v curl > /dev/null 2>&1; then
	curl -fsSL -o "$tmp" "$url" || true
    else
	wget -q -O "$tmp" "$url" || true
    fi

    if [ ! -s "$tmp" ]; then
	rm -f "$tmp"
	echo "$0: could not download '$table' from $url." >&2
	exit 1
    fi

    # A share that serves an HTML page instead of the file, or a file that
    # was replaced, fails here rather than at the first physics result.
    if command -v sha256sum > /dev/null 2>&1; then
	got=$(sha256sum "$tmp" | cut -d' ' -f1)
	if [ "$got" != "$sha256" ]; then
	    rm -f "$tmp"
	    echo "$0: '$table' downloaded from $url has SHA-256 $got," >&2
	    echo "     but $sha256 was expected. Refusing to use it." >&2
	    exit 1
	fi
    else
	echo "$0: sha256sum is not available; the download is unverified."
    fi

    mv "$tmp" "$table"
    exit 0
fi

cat >&2 <<EOF

ERROR: no radiation yields table found for this example.

Looked for:
  ./$table
  $source_file

'$table' is not one of the published tables, so it cannot be downloaded.
Generate it with pychem's pychem_generate_hdf5_parameters, using the
parameter file this example's table was generated from (a piecewise-fits
parameter file gives a mass-only 'M' table; a spectral one gives a mass x
metallicity 'M,Z' table -- match whichever this example's yields_table
name expects), then either copy it here under the name above or point
GEAR_RADIATION_TABLE at it:

  GEAR_RADIATION_TABLE=/path/to/table.hdf5 ./run.sh

The published tables are PopII_parsec_spectral.hdf5, PopIII_parsec_spectral.hdf5
and radiation_fits_popII.hdf5. The tables served by the public hosts (see
getChemistryTable.sh) carry no Data/Radiation group and cannot drive this
example.

EOF
exit 1
