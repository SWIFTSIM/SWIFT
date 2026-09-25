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
#
# For any table name with a published SHA-256 below, steps 1 to 4 all
# verify the file against it, not only the download: a stale or
# hand-replaced local copy is the more likely wrong-table cause and must
# fail loudly too. A table name with no published hash (the
# generate-your-own-with-pychem workflow) is never checked, at any step.
# A mismatch at steps 1 to 3 can be downgraded from a hard failure to a
# warning with GEAR_RADIATION_TABLE_ALLOW_MISMATCH=1, for a deliberately
# hand-edited copy of a published table; the download at step 4 has no
# such override, since a fresh download is never a deliberate edit.

set -eu

table="${1:-PopII_parsec_spectral.hdf5}"

# Published spectral tables: pychem, PARSEC stellar evolution, spectral
# Q_H and L_PE/L_LW. The checksum is what makes a truncated download, a
# silently replaced share, or a stale local file fail loudly instead of
# running.
case "$table" in
    PopII_parsec_spectral.hdf5)
	share_id="ydicQptiff7WspX"
	sha256="b8ba64e393f606b3e373de9d29e5e69c39e5621f3aaa6cd509d1de4cc63e958d"
	;;
    PopIII_parsec_spectral.hdf5)
	share_id="9D5yQEAa3NNg4F8"
	sha256="2de39a002380ba67a8cc931aec59f20c82bab1ec326369d33af5c1fe1d7859d3"
	;;
    # Mass-only ("M"), blackbody-fits Q_H/L table. The HIIRegions and
    # RadiationPressure examples' star masses are calibrated against this
    # table's own Q_H (e.g. Starbench's 26.75 Msun reproduces Bisbas et
    # al. 2015's 1e49 photons/s), so those examples require this exact
    # file, not the spectral one above.
    radiation_fits_popII.hdf5)
	share_id="2bdJyajQjwHnEeB"
	sha256="5594f359bf3431c1832d44aa609767f74349313c35f9cb734d47a8cdcde201a6"
	;;
    *)
	share_id=""
	sha256=""
	;;
esac

# Verify $1 (a file about to become, or already, $table) against the
# published pin for $table, if any. $2 describes where the file came
# from, for the messages. $3 is the specific remedy to print on a
# mismatch. $4 is "hatch" to let GEAR_RADIATION_TABLE_ALLOW_MISMATCH
# downgrade a mismatch to a warning, or "no-hatch" to keep it an
# unconditional hard failure. $5, if given, is a temp file to remove
# before an unconditional exit. Returns normally when the file is
# unpinned, unverifiable (no sha256sum), or matches; exits the whole
# script on an un-escaped mismatch.
verify_table_checksum() {
    local file="$1" origin="$2" fix="$3" hatch="$4" cleanup="${5:-}"

    if [ -z "$sha256" ]; then
	# Not one of the published tables: no pin to check against, this
	# is the documented generate-your-own-with-pychem workflow.
	return 0
    fi

    if ! command -v sha256sum > /dev/null 2>&1; then
	echo "$0: sha256sum is not available; '$table' $origin is unverified." >&2
	return 0
    fi

    local got
    got=$(sha256sum "$file" | cut -d' ' -f1)
    if [ "$got" = "$sha256" ]; then
	return 0
    fi

    echo "$0: '$table' $origin has SHA-256 $got," >&2
    echo "     but $sha256 was expected for the published table of that name." >&2
    echo "     $fix" >&2

    if [ "$hatch" = "hatch" ] && [ "${GEAR_RADIATION_TABLE_ALLOW_MISMATCH:-0}" = "1" ]; then
	echo "     GEAR_RADIATION_TABLE_ALLOW_MISMATCH is set: proceeding anyway." >&2
	return 0
    fi

    echo "     Refusing to use it." >&2
    if [ "$hatch" = "hatch" ]; then
	echo "     If this is a deliberate local edit of a published table," >&2
	echo "     set GEAR_RADIATION_TABLE_ALLOW_MISMATCH=1 to use it anyway." >&2
    fi

    [ -n "$cleanup" ] && rm -f "$cleanup"
    exit 1
}

if [ -e "$table" ]; then
    verify_table_checksum "$table" "already present in this directory" \
	"Remove or rename ./$table and re-run to fetch the published copy." \
	hatch
    exit 0
fi

if [ -n "${GEAR_RADIATION_TABLE:-}" ]; then
    source_file="$GEAR_RADIATION_TABLE"
    source_fix="Point GEAR_RADIATION_TABLE at a matching copy."
else
    source_file="$HOME/programs/pychem/$table"
    source_fix="If this is a newer pychem build, this script's pin needs re-pinning; otherwise point GEAR_RADIATION_TABLE at a matching copy."
fi

if [ -e "$source_file" ]; then
    verify_table_checksum "$source_file" "at '$source_file'" "$source_fix" hatch
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
    # No GEAR_RADIATION_TABLE_ALLOW_MISMATCH override on this path: a
    # fresh download is never a deliberate local edit.
    verify_table_checksum "$tmp" "downloaded from $url" \
	"The share may be serving a truncated file or an HTML error page; retry, or check the share." \
	no-hatch "$tmp"

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
