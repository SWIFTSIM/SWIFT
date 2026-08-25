#!/bin/bash

# make run.sh fail if a subcommand fails. pipefail matters here specifically:
# swift's own exit code is piped into `tee output.log`, and without it a
# crashed/errored swift run still lets the pipeline "succeed" (tee's own
# exit code), silently producing a run_name output directory that looks
# complete but only contains a startup-error log.
set -eo pipefail

n_threads=${n_threads:=8}  #Number of threads to use
gas_density=${gas_density:=1e4} #Gas density in atom/cm^3; see README's Density and metallicity section
gas_particle_mass=${gas_mass:=0.1} #Mass of the gas particles (Msun)
star_mass=${star_mass:=29.7} #Star mass (Msun)
star_type=${star_type:="single_star"}
level=${level:=4} #Resolution level: N = (2**level)**3 gas particles
time_end=${time_end:=1e-4} #TimeIntegration:time_end override (internal units)
dt_max=${dt_max:=1e-5} #TimeIntegration:dt_max override (internal units)
run_name=${run_name:=""}

# Remove the ICs
if [ -e ICs_radiation_pressure.hdf5 ]
then
    rm ICs_radiation_pressure.hdf5
fi

echo "Generating initial conditions to run the example..."
python3 makeIC.py --level $level --rho $gas_density \
	--mass $gas_particle_mass --star_mass $star_mass \
	--star_type $star_type \
    -o ICs_radiation_pressure.hdf5

# Create output directory
DIR=snap
if [ -d "$DIR" ];
then
    echo "$DIR directory exists. Its content will be removed."
    rm -r $DIR
else
    echo "$DIR directory does not exists. It will be created."
    mkdir $DIR
fi

printf "Running simulation..."

# --external-gravity with no Potential: block gives zero gravitational
# force. No --cooling: see README's Configure section.
../../../../swift --hydro --stars --external-gravity --feedback \
		   --sync --limiter --verbose=0 --threads=$n_threads \
		   -P TimeIntegration:time_end:$time_end \
		   -P TimeIntegration:dt_max:$dt_max \
		   params.yml 2>&1 | tee output.log

# Independently recompute the expected momentum output and compare
# against the gas's actual momentum gain
python3 radiation_pressure_momentum_check.py

if [ -z "$run_name" ]; then
    echo "run_name is empty."
else
    if [ -d "$run_name" ]; then
	echo "$run_name directory exists. Nothing will be moved."
    else
	echo "$run_name directory does not exists. It will be created."
	mkdir -p $run_name
	mv snap $run_name
	mv output.log $run_name
	mv timesteps.txt $run_name
	mv statistics.txt $run_name
	mv unused_parameters.yml $run_name
	mv used_parameters.yml $run_name
    fi
fi
