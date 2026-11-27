#!/bin/bash

# Shared GEAR example scripts (tables, plotting)
scripts_location="../../../../GEAR_ICs_and_SCRIPTS"

# swift's exit status is piped into tee: without pipefail a refused or crashed run looks complete.
set -eo pipefail

leg=${leg:="main"}   # main: the cooling reads the kernel mean of the light; own: the particle's own light
n_threads=${n_threads:=8}   # Number of threads
level=${level:=5}   # (2**level)**3 gas particles per glass cube, two cubes stacked in z
gas_density=${gas_density:=1316}   # Gas density in atom/cm^3 (n_H = 1000 cm^-3)
gas_mass=${gas_mass:=0.03825}   # Gas particle mass in Msun: 0.03825 at level 5 gives a 3.38 pc box side
star_mass=${star_mass:=4.9}   # Star mass in Msun
run_name=${run_name:=$leg}
ic_file=ICs_isrf_h2_sheet.hdf5
params=params.yml

case "$leg" in
    main) own=0 ;;
    own) own=1 ;;
    *) echo "Unknown leg '$leg'." ; exit 1 ;;
esac

glass_file="glassCube_$((2**level)).hdf5"
if [ ! -e "$glass_file" ]; then
    echo "Fetching the glass file ($glass_file)..."
    ./getGlass.sh $((2**level))
fi
if [ ! -e CloudyData_UVB=HM2012.h5 ]; then
    "$scripts_location"/getGrackleCoolingTable.sh
fi
"$scripts_location"/getRadiationTable.sh PopII_parsec_spectral.hdf5 || exit 1
"$scripts_location"/checkRadiationTable.sh PopII_parsec_spectral.hdf5 --with-isrf || exit 1

rm -f "$ic_file"
python3 makeIC.py --level $level --rho $gas_density --mass $gas_mass --star_mass $star_mass -o "$ic_file"

rm -rf snap
mkdir snap
swift_exe=${swift_exe:=../../../../../swift}

set +e
$swift_exe --hydro --stars --external-gravity --feedback --cooling \
    --sync --limiter --verbose=0 --threads=$n_threads \
    -P GEARFeedback:ISRF_cooling_reads_own_field:$own \
    $params 2>&1 | tee output.log
swift_status=${PIPESTATUS[0]}
set -e
[ $swift_status -eq 0 ] || { echo "swift stopped with status $swift_status"; exit $swift_status; }

mkdir -p "$run_name"
cp used_parameters.yml output.log timesteps.txt fixture.json "$run_name"/
rm -rf "$run_name"/snap
mv snap "$run_name"/
