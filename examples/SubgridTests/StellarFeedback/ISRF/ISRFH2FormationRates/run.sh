#!/bin/bash

# Shared GEAR example scripts (tables)
scripts_location="../../../../GEAR_ICs_and_SCRIPTS"
starlit=../ISRFH2StarLit

# swift's exit status is piped into tee: without pipefail a refused or crashed run looks complete.
set -eo pipefail

arm=${arm:="dust"}   # hminus: H- channel only (zero metallicity), dust: dust channel only (0.1 solar metallicity)
n_threads=${n_threads:=4}   # Number of threads
level=${level:=4}   # Resolution level: N = (2**level)**3 gas particles
gas_density=${gas_density:=1315.8}   # Gas density in atom/cm^3 (n_H = 0.76 times this, 1000 cm^-3)
gas_mass=${gas_mass:=0.1}   # Gas particle mass in Msun
extra=${extra:=""}   # Extra swift arguments, for example "-P GeajayCooling:reference_dust_to_gas_ratio:0.018774"
run_name=${run_name:=$arm}
force=${force:=0}   # 1: remove an existing snap/ before the run
ic_file=ICs_isrf_h2_formation.hdf5

case "$arm" in
    hminus) metallicity=0; hii_ratio=1e-5 ;;
    dust)   metallicity=0.1; hii_ratio=1e-8 ;;
    *) echo "Unknown arm '$arm'. Use hminus or dust." >&2; exit 1 ;;
esac

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

rm -f "$ic_file"
python3 makeIC.py --level $level --rho $gas_density --mass $gas_mass -o "$ic_file"

rm -rf snap
mkdir snap
../../../../../swift --hydro --stars --external-gravity --feedback --cooling \
    --sync --limiter --verbose=0 --threads=$n_threads \
    -P GEARChemistry:initial_metallicity:$metallicity \
    -P GeajayCooling:initial_hii_to_h_ratio:$hii_ratio $extra params.yml 2>&1 | tee output.log

mkdir -p "$run_name"
cp used_parameters.yml output.log timesteps.txt "$run_name"/
mv snap "$run_name"/

echo "Check the run:"
echo "python3 check_h2_formation_rates.py --arm $arm $run_name/snap"
