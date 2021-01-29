#!/bin/bash

set -e

# MUSIC binary
music=~/music/MUSIC
if test -f $music; then
    echo "Using the following version of MUSIC $music."
else
    echo "MUSIC is not found."
    exit
fi

# Grab the cooling and yield tables if they are not present.
if [ ! -e CloudyData_UVB=HM2012_shielded.h5 ]
then
    echo "Fetching tables..."
    ../getChemistryTable.sh
    ../../Cooling/getGrackleCoolingTable.sh
fi

echo "Generating the initial conditions"
$music music.conf

echo "Converting the initial conditions into a SWIFT compatible format"
python3 convert_ic.py

echo "Running SWIFT"
../../swift --cooling --feedback --cosmology  --limiter --sync --self-gravity --hydro --stars --star-formation --threads=24 agora_cosmo.yml 2>&1 | tee output.log
