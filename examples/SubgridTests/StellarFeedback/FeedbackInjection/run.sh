#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

n_ranks=${n_ranks:=0}  # Number of ranks to use
n_threads=${n_threads:=8}  #Number of threads to use
level=${level:=5}  #Number of particles = 2^(3*level)
gas_density=${gas_density:=0.1} #Gas density in atom/cm^3
gas_particle_mass=${gas_particle_mass:=10} #Mass of the gas particles
scale_height=${z_scale:=0.1} #Mass of the gas particles
star_mass=${star_mass:=29.7} #Mass of the gas particles
seed=${seed:=1} # Random seed for random numbers

ics_filename="ICs_disc.hdf5"

# Remove the ICs
if [ -e $ics_filename ]
then
    rm $ics_filename
fi

#Create the ICs if they do not exist
if [ ! -e $ics_filename ]
then
    echo "Generating initial conditions to run the example..."
    python3 makeIC.py --level $level --rho $gas_density \
	    --mass $gas_particle_mass --star_mass $star_mass \
	    --z_scale $scale_height --seed $seed \
	    -o $ics_filename
fi

# Get the Grackle cooling table
if [ ! -e CloudyData_UVB=HM2012.h5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    ./getGrackleCoolingTable.sh
fi

if [ ! -e POPIIsw.h5 ]
then
    echo "Fetching the chemistry tables..."
    ./getChemistryTable.sh
fi

# Create output directory
DIR="snap"
if [ -d "$DIR" ]; then
    echo "$DIR directory exists. Its content will be removed."
    rm -r "$DIR"
fi
mkdir "$DIR"

if [[ $n_ranks -gt 0 ]]; then
  swift="mpirun -n $n_ranks ../../../../swift_mpi"
else
  swift="../../../../swift"
fi
    
printf "Running simulation..."
$swift --hydro --stars --external-gravity --feedback \
	       --sync --limiter --pin --threads=$n_threads \
	       params.yml 2>&1 | tee output.log
