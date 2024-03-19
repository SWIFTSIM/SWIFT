#!/bin/bash

# exit if anything fails
set -e
set -o pipefail

 # Generate the initial conditions if they are not present.
if [ ! -e glassPlane_128.hdf5 ]
then
    echo "Fetching initial glass file for the 2D RT advection example..."
    ./getGlass.sh
fi
if [ ! -e advection_2D.hdf5 ]
then
    echo "Generating initial conditions for the 2D RT advection example..."
    python3 makeIC.py
fi

# Default run
ymlfile=rt_advection2D_high_redshift.yml
zdomain="h"

# Do we have a cmdline argument provided?
if [ $# -gt 0 ]; then
    case "$1" in
    -l | -low | --l | --low | l | ./rt_advection2D_low_redshift.yml | rt_advection2D_low_redshift | rt_advection2D_low_redshift.yml )
        ymlfile=rt_advection2D_low_redshift.yml
	zdomain="l"
        ;;
    -m | -mid | --m | --mid | m | ./rt_advection2D_medium_redshift.yml | rt_advection2D_medium_redshift | rt_advection2D_medium_redshift.yml )
        ymlfile=rt_advection2D_medium_redshift.yml
	zdomain="m"
        ;;
    -h | -high | -hi | --h | --hi | --high | h | ./rt_advection2D_high_redshift.yml | rt_advection2D_high_redshift | rt_advection2D_high_redshift.yml )
        ymlfile=rt_advection2D_high_redshift.yml
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
