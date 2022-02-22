#!/bin/bash

# make run.sh fail if a subcommand fails
set -e
set -o pipefail

if [ ! -f 'uniformBox-rt.hdf5' ]; then
    echo "Generating ICs"
    python3 makeIC.py
fi

cmd=../../swift
if [ $# -gt 0 ]; then
    case "$1" in 
    g | gdb)
        cmd='gdb --args ../../swift'
        ;;
    *)
        echo unknown cmdline param, running without gdb
        ;;
    esac
fi

# Run SWIFT with RT
$cmd \
    --hydro --threads=4 --stars --external-gravity \
    --feedback --radiation \
    uniform_rt_timestep_output_sync.yml 2>&1 | tee output.log

echo "running sanity checks"
python3 ./rt_sanity_checks.py | tee sanity_check.log
echo "running checks for uniform box test"
python3 ./rt_uniform_box_checks.py | tee box_check.log
