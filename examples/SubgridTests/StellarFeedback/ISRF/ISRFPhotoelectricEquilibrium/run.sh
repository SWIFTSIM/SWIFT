#!/bin/bash

# Shared GEAR example scripts (tables)
scripts_location=${scripts_location:="../../../../GEAR_ICs_and_SCRIPTS"}
swift_bin=${swift_bin:="../../../../../swift"}

# swift's exit status is piped into tee: pipefail keeps a failed run from looking complete.
set -eo pipefail

n_threads=${n_threads:=8}                      # Number of threads to use
gas_density=${gas_density:=10}                 # Gas density in H atoms/cm^3
g0=${g0:=1.65}                                 # Far-UV field in Habing units (1.65 stands for chi = 1 of the Draine field, approximately)
initial_temperature=${initial_temperature:=1000} # Starting temperature in K; use a hot and a cold start to find both branches
gas_particle_mass=${gas_mass:=0.1}             # Mass of the gas particles (Msun)
level=${level:=5}                              # Resolution level: N = (2**level)**3 gas particles
metallicity=${metallicity:=0.6475}             # GEARChemistry:initial_metallicity; 0.6475 is a metal mass fraction of 0.01295
field=${field:=seeded}                         # seeded: field from the ICs; none: no field (control); background: geajay constant background
duration_myr=${duration_myr:=""}               # Run length in Myr; empty: 5 / gas_density
pe_model=${pe_model:=""}                       # Heating model of the cooling library, empty keeps the params file value
run_name=${run_name:=""}
extra_args=${extra_args:=""}                   # Further swift arguments, e.g. "-P Snapshots:delta_time:1e-5"

# The cooling library is read from the binary's configure line.
cooling_library=$(strings "$swift_bin" | grep -o 'with-cooling=[a-z_0-9]*' | head -1 | cut -d= -f2)
params_file=params.yml
pe_key=GrackleCooling:photoelectric_heating_efficiency
case "$cooling_library" in
    geajay_*)
        params_file=params_geajay.yml
        pe_key=GeajayCooling:photoelectric_heating_model
        ;;
    grackle_0)
        echo "The photoelectric heating needs a cooling mode that tracks species:" \
             "grackle_0 has no electron density and the Grackle grain recombination" \
             "term returns NaN for cold gas there. Build with grackle_1 or higher," \
             "or with geajay." >&2
        exit 1
        ;;
esac
pe_arg=""
if [ -n "$pe_model" ]; then
    pe_arg="-P $pe_key:$pe_model"
fi

# Field handling. The ICs always hold the seeded field; "none" writes zero.
ic_g0=$g0
isrf_switch=1
field_args=""
case "$field" in
    seeded) ;;
    none) ic_g0=0 ;;
    background)
        case "$cooling_library" in
            geajay_*) ;;
            *) echo "field=background needs a geajay build." >&2; exit 1 ;;
        esac
        ic_g0=0
        isrf_switch=0
        field_args="-P GeajayCooling:background_mode:per_band \
                    -P GeajayCooling:ionising_background_model:none \
                    -P GeajayCooling:lw_background_model:none \
                    -P GeajayCooling:isrf_background_model:constant \
                    -P GeajayCooling:isrf_background_habing:$g0 \
                    -P GeajayCooling:dust_model:basic"
        ;;
    *) echo "field must be seeded, none or background." >&2; exit 1 ;;
esac

# Time step and run length from the density. The relaxation time scales as 1/n:
# 0.0064 Myr at 10 H/cm^3 for the constant efficiency. dt_max stays below a sixth of it.
if [ -z "$duration_myr" ]; then
    duration_myr=$(python3 -c "print(5/$gas_density)")
fi
time_end=$(python3 -c "print($duration_myr*1.02271e-3)")
dt_max=$(python3 -c "print(1e-5/$gas_density)")
delta_time=$(python3 -c "print($time_end/20)")
stat_time=$(python3 -c "print($time_end/200)")

[ -e ICs_isrf_photoelectric_equilibrium.hdf5 ] && rm ICs_isrf_photoelectric_equilibrium.hdf5

glass_n=$((2**level))
glass_file="glassCube_${glass_n}.hdf5"
if [ ! -e "$glass_file" ]; then
    echo "Fetching initial glass file (${glass_file})..."
    ./getGlass.sh $glass_n
fi

if [ ! -e CloudyData_UVB=HM2012.h5 ]; then
    echo "Fetching the Cloudy tables required by Grackle..."
    "$scripts_location"/getGrackleCoolingTable.sh
fi

"$scripts_location"/getRadiationTable.sh PopII_parsec_spectral.hdf5 || exit 1

echo "Generating initial conditions to run the example..."
python3 makeIC.py --level $level --rho $gas_density --mass $gas_particle_mass \
    --g0 $ic_g0 -o ICs_isrf_photoelectric_equilibrium.hdf5

DIR=snap
if [ -d "$DIR" ]; then
    echo "$DIR directory exists. Its content will be removed."
    rm -r $DIR
fi
mkdir $DIR

printf "Running simulation..."

# --external-gravity with no Potential: block gives zero force. The field of the ICs is
# never rewritten because there is no star and GEARFeedback:ISRF_propagation is 0.
"$swift_bin" --hydro --stars --external-gravity --feedback --cooling \
    --sync --limiter --verbose=0 --threads=$n_threads \
    -P TimeIntegration:time_end:$time_end \
    -P TimeIntegration:dt_max:$dt_max \
    -P Snapshots:delta_time:$delta_time \
    -P Statistics:delta_time:$stat_time \
    -P SPH:initial_temperature:$initial_temperature \
    -P GEARChemistry:initial_metallicity:$metallicity \
    -P GEARFeedback:with_interstellar_radiation_field:$isrf_switch \
    $field_args $pe_arg $extra_args $params_file 2>&1 | tee output.log

if [ -n "$run_name" ]; then
    if [ -d "$run_name" ]; then
        echo "$run_name directory exists. Nothing will be moved."
    else
        mkdir -p $run_name
        mv snap output.log timesteps.txt statistics.txt unused_parameters.yml used_parameters.yml $run_name
    fi
fi
