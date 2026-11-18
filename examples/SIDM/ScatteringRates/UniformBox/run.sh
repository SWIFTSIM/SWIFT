#!/bin/bash

# SIDM scattering rate in a uniform-density periodic box (cosmological run).
#
# Usage: ./run.sh [swift binary]   (default ../../../../swift)
#
# The same ICs (comoving positions, peculiar velocities) are started at each
# scale factor in A_VALUES. The rates are computed during start-up and
# written to snapshot_0000, so one step per run is enough. Gravity is off.
# check_uniform_box.py then checks that
#
#   1. the rate scales as a^-3 and matches the analytic rate at every a, and
#   2. the kernel densities and rates are biased high by the expected amount
#      once the self-contribution is removed.
#
# Set ETA to change BasicSIDM:resolution_eta. SWIFT must be configured with
# --with-sidm=Basic and GSL (needed for the cosmology tables).

set -e

swift=${1:-../../../../swift}
A_VALUES=${A_VALUES:-"1.0 0.5 0.25 0.1 0.05"}

if [ ! -e uniform_box.hdf5 ]; then
    echo "Generating initial conditions..."
    python3 makeIC.py
fi

for a in $A_VALUES; do
    a_end=$(python3 -c "print($a * 1.0002)")
    printf "Running a = %s...\n" "$a"
    "$swift" --sidm --cosmology -n 1 --threads=8 \
        -P Cosmology:a_begin:$a -P Cosmology:a_end:$a_end \
        -P Snapshots:scale_factor_first:$a -P Snapshots:subdir:snap_a$a \
        -P Statistics:scale_factor_first:$a \
        ${ETA:+-P BasicSIDM:resolution_eta:$ETA} \
        params.yml > output_a$a.log 2>&1
done

printf "Checking scattering rates...\n"
python3 check_uniform_box.py snap_a*/snapshot_0000.hdf5
