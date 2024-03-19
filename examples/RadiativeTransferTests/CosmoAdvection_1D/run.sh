#!/bin/bash

# make run.sh fail if a subcommand fails
set -e
set -o pipefail

if [ ! -f advection_1D.hdf5 ]; then
    echo "Generating ICs"
    python3 makeIC.py
fi

# Default run
ymlfile=rt_advection1D_high_redshift.yml
zdomain="h"

# Do we have a cmdline argument provided?
if [ $# -gt 0 ]; then
    case "$1" in 
    -l | -low | --l | --low | l | ./rt_advection1D_low_redshift.yml | rt_advection1D_low_redshift | rt_advection1D_low_redshift.yml )
        ymlfile=rt_advection1D_low_redshift.yml
	zdomain="l"
        ;;
    -m | -mid | --m | --mid | m | ./rt_advection1D_medium_redshift.yml | rt_advection1D_medium_redshift | rt_advection1D_medium_redshift.yml )
        ymlfile=rt_advection1D_medium_redshift.yml
	zdomain="m"
        ;;
    -h | -high | -hi | --h | --hi | --high | h | ./rt_advection1D_high_redshift.yml | rt_advection1D_high_redshift | rt_advection1D_high_redshift.yml )
        ymlfile=rt_advection1D_high_redshift.yml
	zdomain="h"
        ;;
    *)
        echo unknown cmdline param, running default $ymlfile
        ;;
    esac
fi



# Run SWIFT with RT
../../../swift \
    --hydro \
    --cosmology \
    --threads=4 \
    --verbose=0  \
    --radiation \
    --stars \
    --feedback \
    --external-gravity \
    $ymlfile 2>&1 | tee output.log

python3 ./plotSolution.py -z $zdomain
python3 ./plotEnergy.py -z $zdomain
