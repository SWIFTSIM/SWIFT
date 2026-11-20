#!/bin/bash

# Shared GEAR example scripts (tables)
scripts_location="../../../../GEAR_ICs_and_SCRIPTS"
starlit=../ISRFH2StarLit

# swift's exit status is piped into tee: without pipefail a refused or crashed run looks complete.
set -eo pipefail

leg=${leg:="on"}   # on: line absorption on, off: line absorption off, repeat: a second run of the same switch
switch=${switch:=$([ "$leg" = "off" ] && echo 0 || echo 1)}
n_threads=${n_threads:=8}   # Number of threads
level=${level:=5}   # Resolution level: N = (2**level)**3 gas particles
gas_density=${gas_density:=39.5}   # Gas density in atom/cm^3 (n_H = 0.76 times this)
gas_mass=${gas_mass:=0.8}   # Gas particle mass in Msun (box side about 30 pc)
h2_ratio=${h2_ratio:=1e-4}   # Initial nH2 / nH
run_name=${run_name:=$leg}
ic_file=ICs_isrf_h2_cooling.hdf5

glass_file="glassCube_$((2**level)).hdf5"
if [ ! -e "$glass_file" ]; then
    echo "Fetching the glass file ($glass_file)..."
    "$starlit"/getGlass.sh $((2**level))
fi
if [ ! -e CloudyData_UVB=HM2012.h5 ]; then
    "$scripts_location"/getGrackleCoolingTable.sh
fi
"$scripts_location"/getRadiationTable.sh PopII_parsec_spectral.hdf5 || exit 1
"$scripts_location"/checkRadiationTable.sh PopII_parsec_spectral.hdf5 --with-isrf || exit 1

rm -f "$ic_file"
python3 "$starlit"/makeIC.py --level $level --rho $gas_density --mass $gas_mass \
    --h2_ratio $h2_ratio -o "$ic_file"

rm -rf snap
mkdir snap
../../../../../swift --hydro --stars --external-gravity --feedback --cooling \
    --sync --limiter --verbose=0 --threads=$n_threads \
    -P GEARFeedback:ISRF_H2_line_absorption:$switch params.yml 2>&1 | tee output.log

mkdir -p "$run_name"
cp used_parameters.yml output.log timesteps.txt fixture.json "$run_name"/
mv snap "$run_name"/

echo "Run the legs on and off, then compare them:"
echo "python3 check_h2_cooling_onoff.py --on on/snap --off off/snap --first 3 --last 6 --noise <value>"
