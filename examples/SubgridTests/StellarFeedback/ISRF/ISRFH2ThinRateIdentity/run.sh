#!/bin/bash

# Shared GEAR example scripts (tables)
scripts_location="../../../../GEAR_ICs_and_SCRIPTS"
starlit=../ISRFH2StarLit

# swift's exit status is piped into tee: without pipefail a refused or crashed run looks complete.
set -eo pipefail

n_threads=${n_threads:=4}   # Number of threads
level=${level:=5}   # Resolution level: N = (2**level)**3 gas particles
gas_density=${gas_density:=1e3}   # Gas density in atom/cm^3
gas_mass=${gas_mass:=0.1}   # Gas particle mass in Msun
star_mass=${star_mass:=29.7}   # Star mass in Msun
h2_ratio=1e-6   # Initial nH2 / nH; fixed, it must equal GeajayCooling:initial_h2i_to_h_ratio of params.yml
extra=${extra:=""}   # Extra swift arguments
run_name=${run_name:=thin}
force=${force:=0}   # 1: remove an existing snap/ and $run_name/ before the run
ic_file=ICs_isrf_h2_thin_absolute.hdf5

if [ "$force" != 1 ] && { [ -e snap ] || [ -e "$run_name" ]; }; then
    echo "snap/ or $run_name/ exists: move it away, or run with force=1 to remove it." >&2
    exit 1
fi
rm -rf "$run_name"

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
    --star_mass $star_mass --h2_ratio $h2_ratio -o "$ic_file"

rm -rf snap
mkdir snap
../../../../../swift --hydro --stars --external-gravity --feedback --cooling \
    --sync --limiter --verbose=0 --threads=$n_threads $extra params.yml 2>&1 | tee output.log

mkdir -p "$run_name"
cp used_parameters.yml output.log timesteps.txt fixture.json "$run_name"/
mv snap "$run_name"/

echo "Check the run:"
echo "python3 check_h2_thin_absolute.py $run_name/snap --radiation-table PopII_parsec_spectral.hdf5"
