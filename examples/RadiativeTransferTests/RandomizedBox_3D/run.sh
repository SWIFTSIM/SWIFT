#!/bin/bash

# make run.sh fail if a subcommand fails
set -e
set -o pipefail

if [ ! -f 'randomized-sine.hdf5' ]; then
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
    --hydro \
    --threads=9 \
    --verbose=0  \
    --radiation \
    --self-gravity \
    --stars \
    --feedback \
    ./randomized-rt.yml 2>&1 | tee output.log

echo "running sanity checks"
python3 ../UniformBox_3D/rt_sanity_checks.py | tee sanity_check.log
