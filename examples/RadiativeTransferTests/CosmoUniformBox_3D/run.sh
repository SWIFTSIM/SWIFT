#!/bin/bash

# exit if anything fails
set -e
set -o pipefail

# Generate the initial conditions if they are not present.
if [ ! -e uniform_3D.hdf5 ]
then
    echo "Generating initial conditions for the 3D uniform box example..."
    python3 makeIC.py
fi

# Default run
ymlfile=rt_uniform3D_high_redshift.yml
zdomain="h"

# Do we have a cmdline argument provided?
if [ $# -gt 0 ]; then
    case "$1" in
    -l | -low | --l | --low | l | ./rt_uniform3D_low_redshift.yml | rt_uniform3D_low_redshift | rt_uniform3D_low_redshift.yml )
        ymlfile=rt_uniform3D_low_redshift.yml
	zdomain="l"
        ;;
    -m | -mid | --m | --mid | m | ./rt_uniform3D_medium_redshift.yml | rt_uniform3D_medium_redshift | rt_uniform3D_medium_redshift.yml )
        ymlfile=rt_uniform3D_medium_redshift.yml
	zdomain="m"
        ;;
    -h | -high | -hi | --h | --hi | --high | h | ./rt_uniform3D_high_redshift.yml | rt_uniform3D_high_redshift | rt_uniform3D_high_redshift.yml )
        ymlfile=rt_uniform3D_high_redshift.yml
	zdomain="h"
        ;;
    *)
        echo unknown cmdline param, running default $ymlfile
        ;;
    esac
fi


# Run SWIFT with RT and cosmology
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
