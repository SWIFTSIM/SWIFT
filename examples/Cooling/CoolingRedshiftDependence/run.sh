#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e gravity_glassCube_32.hdf5 ]
then
    echo "Fetching initial gravity glass file for the constant cosmological box example..."
    ./getGlass.sh
fi

# Fetch the cooling tables
if [ ! -e ics_no_z.hdf5 ]
then
    echo "Generating initial conditions for the uniform cosmo box example..."
    python3 makeIC.py
fi

swift_location="../../../swift"

rm data/redshift_dependence_*_z_*.hdf5

mkdir data

# Run SWIFT
$swift_location --hydro --cosmology --cooling --limiter --threads=4 cooling_redshift_dependence_high_z.yml 2>&1 | tee output_high.log
mv timesteps.txt timesteps_high.txt
$swift_location --hydro --cosmology --cooling --limiter --threads=4 cooling_redshift_dependence_low_z.yml 2>&1 | tee output_low.log
mv timesteps.txt timesteps_low.txt
$swift_location --hydro --cooling --limiter --threads=4 cooling_redshift_dependence_no_z.yml 2>&1 | tee output_no.log
mv timesteps.txt timesteps_no.txt

python3 plotSolution.py
