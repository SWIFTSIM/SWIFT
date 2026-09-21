#!/bin/bash
#
# Check that a GEAR yields table can drive the radiation feedback.
#
# Usage: checkRadiationTable.sh <table.h5> [--with-isrf]
#
# Any example running with photoionization, radiation pressure or the
# interstellar radiation field needs a table carrying a "Data/Radiation"
# group. With --with-isrf it must also carry the four ISRF band datasets.
# The tables served by the public hosts (see getChemistryTable.sh) predate
# the radiation tables and carry neither, so this check stops the example
# before it spends time on the glass file, the Cloudy tables and the
# initial conditions.

set -u

table="${1:-}"
with_isrf=0

if [ -z "$table" ]; then
    echo "Usage: $0 <table.h5> [--with-isrf]" >&2
    exit 2
fi
shift

while [ $# -gt 0 ]; do
    case "$1" in
	--with-isrf)
	    with_isrf=1
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

python3 - "$table" "$with_isrf" <<'EOF'
import sys

import h5py

table, with_isrf = sys.argv[1], sys.argv[2] == "1"
isrf_fields = ("L_FUV", "L_LW", "Integrated_L_FUV", "Integrated_L_LW")

with h5py.File(table, "r") as f:
    if "Data/Radiation" not in f:
        missing = "the 'Data/Radiation' group"
    else:
        absent = [d for d in isrf_fields if d not in f["Data/Radiation"]]
        missing = None
        if with_isrf and absent:
            missing = "the dataset(s) " + ", ".join(
                "'Data/Radiation/%s'" % d for d in absent
            )

if missing is None:
    sys.exit(0)

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
sys.exit(1)
EOF
