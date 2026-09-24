#!/bin/bash
#
# Check that a GEAR yields table can drive the radiation feedback.
#
# Usage: checkRadiationTable.sh <table.h5> [--with-isrf] [--require-1d]
#
# Any example running with photoionization, radiation pressure or the
# interstellar radiation field needs a table carrying a "Data/Radiation"
# group, and that group must carry "MeanPhotonEnergyLW"/
# "Integrated_MeanPhotonEnergyLW": SWIFT requires both unconditionally,
# not only under with_interstellar_radiation_field. With --with-isrf it
# must also carry the four ISRF band datasets.
# The tables served by the public hosts (see getChemistryTable.sh) predate
# the radiation tables and carry neither, so this check stops the example
# before it spends time on the glass file, the Cloudy tables and the
# initial conditions.
#
# --require-1d additionally fails on a mass x metallicity ("M,Z") table.
# SWIFT itself reads either dimensionality, but the example-side Python
# analysis (radiation_table_reader.py, used by the Starbench, Hu2017,
# StromgrenSphere(Cosmo) and RadiationPressureShellExpansion check
# scripts) only understands a mass-only ("M") table and raises
# NotImplementedError on the other. Pass this flag from an example whose
# own check script imports that module.
#
# SWIFT itself logs the table's identity when it reads the table, so a
# normal run's log already records which table drove it. This check prints
# the same attributes for the case SWIFT never reaches: a table that fails
# below, where the identity is what says why.

set -u

table="${1:-}"
with_isrf=0
require_1d=0

if [ -z "$table" ]; then
    echo "Usage: $0 <table.h5> [--with-isrf] [--require-1d]" >&2
    exit 2
fi
shift

while [ $# -gt 0 ]; do
    case "$1" in
	--with-isrf)
	    with_isrf=1
	    shift
	    ;;
	--require-1d)
	    require_1d=1
	    shift
	    ;;
	*)
	    echo "$0: unknown option '$1'." >&2
	    exit 2
	    ;;
    esac
done

if [ ! -e "$table" ]; then
    echo "$0: '$table' does not exist." >&2
    exit 1
fi

if ! python3 -c "import h5py" > /dev/null 2>&1; then
    echo "$0: python3 with h5py is not available; skipping the check of" \
	 "'$table'."
    exit 0
fi

python3 - "$table" "$with_isrf" "$require_1d" <<'EOF'
import sys

import h5py

table = sys.argv[1]
with_isrf = sys.argv[2] == "1"
require_1d = sys.argv[3] == "1"
isrf_fields = ("L_PE", "L_LW", "Integrated_L_PE", "Integrated_L_LW")
required_fields = ("MeanPhotonEnergyLW", "Integrated_MeanPhotonEnergyLW")

# Attributes that identify the table itself. GEARFeedback:yields_table only
# names a file, and the name says nothing about which Q_H the run used, so
# print these to the run log. qh_source and stellar_evolution_source are
# carried by the current pychem tables; older tables carry a single source
# attribute instead. mass_min_msun, mass_max_msun and nz separate a Pop II
# table from a Pop III one, which the source attributes alone do not.
identity_keys = (
    "qh_source",
    "lwpe_source",
    "stellar_evolution_source",
    "source",
    "dimensionality",
    "mass_min_msun",
    "mass_max_msun",
    "nz",
)


def as_text(value):
    """Return an HDF5 attribute as printable text, decoding bytes."""
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    return str(value)


identity = None
dimensionality = None

# A file that is not HDF5 at all reaches here: getRadiationTable.sh copies
# whatever GEAR_RADIATION_TABLE names, and takes an existing file in the
# example directory as given. Report it in this script's own words, so the
# run log says which check stopped the example.
try:
    handle = h5py.File(table, "r")
except OSError:
    sys.stderr.write(
        "\n"
        "ERROR: '%s' is not a readable HDF5 file.\n"
        "\n"
        "Stage a GEAR yields table carrying a 'Data/Radiation' group under\n"
        "that name, or point GEAR_RADIATION_TABLE at one.\n"
        "\n" % table
    )
    sys.exit(1)

with handle as f:
    if "Data/Radiation" not in f:
        missing = "the 'Data/Radiation' group"
    else:
        group = f["Data/Radiation"]
        identity = [
            (key, as_text(group.attrs[key]))
            for key in identity_keys
            if key in group.attrs
        ]
        if "dimensionality" in group.attrs:
            dimensionality = as_text(group.attrs["dimensionality"])
        absent = [d for d in required_fields if d not in group]
        if with_isrf:
            absent += [d for d in isrf_fields if d not in group]
        missing = None
        if absent:
            missing = "the dataset(s) " + ", ".join(
                "'Data/Radiation/%s'" % d for d in absent
            )

# Printed before the verdict, so a table that fails the check below still
# says what it is. A table with no identity attribute is reported as such
# and still passes: the provenance is metadata, not something the radiation
# reader needs.
if identity is not None:
    if identity:
        sys.stdout.write("Radiation table '%s' reports:\n" % table)
        for key, value in identity:
            sys.stdout.write("  %s = %s\n" % (key, value))
    else:
        sys.stdout.write(
            "Radiation table '%s' carries no provenance attribute in its "
            "'Data/Radiation' group.\n" % table
        )
    sys.stdout.flush()

# Independent checks: a missing dataset says nothing about dimensionality
# and vice versa, so both are reported when both apply.
dimensionality_fails = require_1d and dimensionality not in (None, "M")
failed = False

if missing is not None:
    failed = True
    sys.stderr.write(
        "\n"
        "ERROR: '%s' is missing %s.\n"
        "\n"
        "This example runs GEAR radiation feedback, which reads the stellar\n"
        "photon rates and luminosities from that group. SWIFT aborts at\n"
        "start-up on a table without it.\n"
        "\n"
        "The tables on the public hosts predate the GEAR radiation tables and\n"
        "do not carry this data. Generate a table with pychem's\n"
        "pychem_generate_hdf5_parameters, or ask the GEAR maintainers for one,\n"
        "then point GEARFeedback:yields_table in params.yml at it.\n"
        "\n" % (table, missing)
    )

if dimensionality_fails:
    failed = True
    sys.stderr.write(
        "\n"
        "ERROR: '%s' is a '%s' (mass x metallicity) table.\n"
        "\n"
        "This example's own check script imports radiation_table_reader.py,\n"
        "which only understands a mass-only ('M') table and raises\n"
        "NotImplementedError on any other dimensionality. SWIFT itself reads\n"
        "either table fine; this is a check-script limitation, not a SWIFT\n"
        "one. Point GEARFeedback:yields_table at a mass-only table instead.\n"
        "\n" % (table, dimensionality)
    )

sys.exit(1 if failed else 0)
EOF
