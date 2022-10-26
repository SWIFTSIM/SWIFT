#!/bin/bash

# make run.sh fail if a subcommand fails
set -e
set -o pipefail

if [ ! -f 'randomized-sine.hdf5' ]; then
    echo "Generating ICs"
    python3 makeIC.py
fi

# use cmdline args as shortcut to run with debugger/MPI
# -g run with gdb
# -m run with MPI
# -mg run with MPI and gdb in individual xterm windows
# -ml run with MPI and write individual output file per MPI rank
cmd=../../../swift
if [ $# -gt 0 ]; then
    case "$1" in 
    -g | g | gdb)
        cmd='gdb --args ../../../swift'
        ;;
    -m | m | mpi)
        cmd='mpirun -n 3 ../../../swift_mpi' 
        ;;
    -mg | -gm | gm | mg | gmpi | gdbmpi )
        cmd='mpirun -n 3 xterm -e gdb -ex run --args ../../../swift_mpi'
        ;;
    -ml | ml | lm | mpilog | logmpi)
        cmd='mpirun -n 2 --output-filename individual_rank_output --merge-stderr-to-stdout ../../../swift_mpi'
        ;;
    *)
        echo unknown cmdline param, running without gdb
        ;;
    esac
fi

# Run SWIFT with RT
$cmd \
    --hydro \
    --threads=3 \
    --verbose=0 \
    --radiation \
    --self-gravity \
    --stars \
    --feedback \
    ./randomized-rt.yml 2>&1 | tee output.log

echo "running sanity checks"
python3 ../UniformBox_3D/rt_sanity_checks.py | tee sanity_check.log
