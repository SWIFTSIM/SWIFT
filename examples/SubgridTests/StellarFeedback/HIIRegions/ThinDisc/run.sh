#!/bin/bash

# Shared GEAR example scripts (tables)
scripts_location="../../../../GEAR_ICs_and_SCRIPTS"
set -eo pipefail

n_threads=${n_threads:=4}           # Number of threads to use
gas_mass=${gas_mass:=20.0}          # Gas particle mass in M_sun
disc_mass=${disc_mass:=4e6}         # Total gas mass in M_sun
scale_length_pc=${scale_length_pc:=150.0} # Exponential radial scale length in parsec
scale_height_pc=${scale_height_pc:=20.0}  # Gaussian vertical sigma in parsec
truncation=${truncation:=5.0}       # Disc radius in scale lengths
boxsize=${boxsize:=1.6}             # Box side in kpc; must exceed 2 * truncation * scale length
n_clusters=${n_clusters:=3}         # Number of star clusters
stars_per_cluster=${stars_per_cluster:=10}
n_field_stars=${n_field_stars:=10}  # Isolated stars
nside=${nside:=1}                   # GEARFeedback:HII_angular_nside override: 0 (spherical) or 1 (12 cones)
random_frame=${random_frame:=0}     # GEARFeedback:HII_random_cone_frame_per_star override:
                                    # 0 = cones fixed to the box axes, 1 = random rotation per star
rebuild_time_myr=${rebuild_time_myr:=0.02} # GEARFeedback:HII_rebuild_time_Myr override
initial_metallicity=${initial_metallicity:=0.02} # GEARChemistry:initial_metallicity override (Z/Zsun)
time_end=${time_end:=3e-3}          # TimeIntegration:time_end override (internal units, 1e-3 = 0.98 Myr)
run_name=${run_name:=""}
restart=${restart:=0}
run_analysis=${run_analysis:=1}     # Run plot_hii_disc.py on the final snapshot

if [ ! -e CloudyData_UVB=HM2012.h5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    "$scripts_location"/getGrackleCoolingTable.sh
fi

"$scripts_location"/getRadiationTable.sh radiation_fits_popII.hdf5 || exit 1
"$scripts_location"/checkRadiationTable.sh radiation_fits_popII.hdf5 || exit 1

rm -f ICs_thin_disc.hdf5
echo "Generating initial conditions to run the example..."
python3 makeIC_disc.py --gas_mass $gas_mass --disc_mass $disc_mass \
        --scale_length_pc $scale_length_pc --scale_height_pc $scale_height_pc \
        --n_clusters $n_clusters --stars_per_cluster $stars_per_cluster \
        --n_field_stars $n_field_stars --truncation $truncation --boxsize $boxsize -o ICs_thin_disc.hdf5

DIR=snap
if [ -d "$DIR" ]; then
    if [ "$restart" -ne 1 ]; then
        echo "$DIR directory exists. Its content will be removed."
        rm -r $DIR
        mkdir $DIR
    fi
else
    mkdir $DIR
fi

if [ "$restart" -eq 1 ]; then
    runtime_param="--restart --verbose=0"
else
    runtime_param="--verbose=0"
fi

echo "Running simulation..."
../../../../../swift --hydro --stars --external-gravity --feedback --cooling \
    --sync --limiter $runtime_param --threads=$n_threads \
    -P GEARFeedback:HII_angular_nside:$nside \
    -P GEARFeedback:HII_random_cone_frame_per_star:$random_frame \
    -P GEARFeedback:HII_rebuild_time_Myr:$rebuild_time_myr \
    -P GEARChemistry:initial_metallicity:$initial_metallicity \
    -P TimeIntegration:time_end:$time_end \
    params.yml 2>&1 | tee output.log

output_dir="."
if [ -n "$run_name" ]; then
    if [ -d "$run_name" ]; then
        echo "$run_name directory exists. Nothing will be moved."
    else
        mkdir -p "$run_name"
        mv snap output.log timesteps.txt statistics.txt unused_parameters.yml used_parameters.yml "$run_name"
        output_dir="$run_name"
    fi
fi

if [ "$run_analysis" -eq 1 ]; then
    echo "Running plot_hii_disc.py on the final snapshot..."
    python3 plot_hii_disc.py -s "$output_dir/snap/snapshot"
fi
