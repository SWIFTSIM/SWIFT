#!/bin/bash

# Shared GEAR example scripts (tables, plotting)
scripts_location=${scripts_location:="../../../../GEAR_ICs_and_SCRIPTS"}
swift_bin=${swift_bin:="../../../../../swift"}

# make run.sh fail if a subcommand fails. pipefail matters here specifically:
# swift's own exit code is piped into `tee output.log`, and without it a
# crashed/errored swift run still lets the pipeline "succeed" (tee's own
# exit code), silently producing a run_name output directory that looks
# complete but only contains a startup-error log.
set -eo pipefail

n_threads=${n_threads:=8}  #Number of threads to use
gas_density=${gas_density:=1e3} #Gas density in atom/cm^3
gas_particle_mass=${gas_mass:=0.1} #Mass of the gas particles (Msun)
star_mass=${star_mass:=29.7} #Star mass (Msun)
star_type=${star_type:="single_star"}
level=${level:=5} #Resolution level: N = (2**level)**3 gas particles
time_end=${time_end:=1.5e-4} #TimeIntegration:time_end override (internal units)
dt_max=${dt_max:=1e-6} #TimeIntegration:dt_max override (internal units)
delta_time=${delta_time:=1e-5} #Snapshots:delta_time override (internal units)
initial_metallicity=${initial_metallicity:=1} #GEARChemistry:initial_metallicity override (Z/Zsun)
with_pe_heating=${with_pe_heating:=1} #GEARFeedback:with_interstellar_radiation_field override (1=on, 0=off)
run_name=${run_name:=""}
pe_model=${pe_model:=""} #Photoelectric heating model of the cooling library (see README); empty keeps the params file value
extra_args=${extra_args:=""} #Further swift arguments, e.g. "-P GEARFeedback:ISRF_propagation:0"

# The cooling library is read from the binary's configure line.
cooling_library=$(strings "$swift_bin" | grep -o 'with-cooling=[a-z_0-9]*' | head -1 | cut -d= -f2 || true)
if [ -z "$cooling_library" ] && strings "$swift_bin" | grep -q 'with-subgrid=GEAR'; then
    cooling_library=grackle_0  # --with-subgrid=GEAR selects grackle_0
fi
if [ -z "$cooling_library" ]; then
    echo "Cannot read the cooling library from $swift_bin: its configure line" \
         "has no with-cooling= or with-subgrid=GEAR. Set swift_bin to a swift" \
         "built with --with-cooling=grackle_N or geajay_N." >&2
    exit 1
fi
params_file=params.yml
pe_key=GrackleCooling:photoelectric_heating_efficiency
case "$cooling_library" in
    geajay_*)
        params_file=params_geajay.yml
        pe_key=GeajayCooling:photoelectric_heating_model
        extra_args="-P GeajayCooling:primordial_chemistry_network:${cooling_library#geajay_} $extra_args"
        ;;
    grackle_0)
        if [ "$with_pe_heating" != "0" ]; then
            echo "Photoelectric heating at grackle_0 gives NaN: the Grackle grain" \
                 "recombination term divides by the electron density, which" \
                 "this mode estimates as zero for cold gas." \
                 "Build with grackle_1 or higher, or with geajay." >&2
            exit 1
        fi
        ;;
esac
pe_arg=""
if [ -n "$pe_model" ]; then
    pe_arg="-P $pe_key:$pe_model"
fi

# Remove the ICs
if [ -e ICs_isrf_photoelectric_heating.hdf5 ]
then
    rm ICs_isrf_photoelectric_heating.hdf5
fi

glass_n=$((2**level))
glass_file="glassCube_${glass_n}.hdf5"
if [ ! -e "$glass_file" ]
then
    echo "Fetching initial glass file (${glass_file})..."
    ./getGlass.sh $glass_n
fi

# Get the Grackle cooling table
if [ ! -e CloudyData_UVB=HM2012.h5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    "$scripts_location"/getGrackleCoolingTable.sh
fi

"$scripts_location"/getRadiationTable.sh PopII_parsec_spectral.hdf5 || exit 1

# Stop here on a table the radiation reader cannot use, rather than
# aborting at start-up once the initial conditions are built.
# The ISRF band datasets are only needed when the run enables the ISRF.
isrf_table_check_flags=""
if [ "$with_pe_heating" != "0" ]; then
    isrf_table_check_flags="--with-isrf"
fi
"$scripts_location"/checkRadiationTable.sh PopII_parsec_spectral.hdf5 \
    $isrf_table_check_flags || exit 1

echo "Generating initial conditions to run the example..."
ic_output=$(python3 makeIC.py --level $level --rho $gas_density \
	--mass $gas_particle_mass --star_mass $star_mass \
	--star_type $star_type \
    -o ICs_isrf_photoelectric_heating.hdf5)
echo "$ic_output"

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
# force. --cooling with GrackleCooling:disable_cooling_for_debugging left
# at 0 (params.yml) means gas actually cools and is heated by Grackle
# every step: see README.
"$swift_bin" --hydro --stars --external-gravity --feedback --cooling \
		   --sync --limiter --verbose=0 --threads=$n_threads \
		   -P TimeIntegration:time_end:$time_end \
		   -P TimeIntegration:dt_max:$dt_max \
		   -P Snapshots:delta_time:$delta_time \
		   -P GEARChemistry:initial_metallicity:$initial_metallicity \
		   -P GEARFeedback:with_interstellar_radiation_field:$with_pe_heating \
		   $pe_arg $extra_args $params_file 2>&1 | tee output.log

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
