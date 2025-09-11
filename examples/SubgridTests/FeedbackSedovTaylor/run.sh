#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

n_threads=${n_threads:=8}  #Number of threads to use
level=${level:=6}  #Number of particles = 2^(3*level)
gas_density=${gas_density:=1} #Gas density in atom/cm^3
gas_particle_mass=${gas_particle_mass:=10} #Mass of the gas particles
star_mass=${star_mass:=29.7} #Mass of the gas particles
with_cooling=${with_cooling:=0}
L=${boxsize:=1} #boxsize in kpc
run_name=${run_name:=""}

# Remove the ICs
if [ -e ICs_homogeneous_box.hdf5 ]
then
    rm ICs_homogeneous_box.hdf5
fi

#Create the ICs if they do not exist
if [ ! -e ICs_homogeneous_box.hdf5 ]
then
    echo "Generating initial conditions to run the example..."
    python3 makeIC.py --level $level --rho $gas_density --boxsize $L \
	    --mass $gas_particle_mass --star_mass $star_mass \
	    -o ICs_homogeneous_box.hdf5
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
DIR=snap #First test of units conversion
if [ -d "$DIR" ];
then
    echo "$DIR directory exists. Its content will be removed."
    rm -r $DIR
else
    echo "$DIR directory does not exists. It will be created."
    mkdir $DIR
fi

printf "Running simulation..."

if [ "$with_cooling" -eq 1 ]; then
../../../swift --hydro --stars --external-gravity --feedback --cooling \
	       --sync --limiter --threads=$n_threads \
	       params.yml 2>&1 | tee output.log
else
../../../swift --hydro --stars --external-gravity --feedback \
	       --sync --limiter --threads=$n_threads \
	       params.yml 2>&1 | tee output.log
fi

#Do some data analysis to show what's in this box
echo "Plotting gas density"
python3 plot_gas_density.py -i 10 -s 'snap/snapshot'
python3 plot_gas_density.py -i 20 -s 'snap/snapshot'
python3 plot_gas_density.py -i 30 -s 'snap/snapshot'
python3 plot_gas_density.py -i 40 -s 'snap/snapshot'
python3 plot_gas_density.py -i 50 -s 'snap/snapshot'
python3 plot_gas_density.py -i 100 -s 'snap/snapshot'

# Movie time
python3 plot_gas_density.py -i 0 -f 100 -s 'snap/snapshot'

# Don't do a phase-space diagram if we are running without cooling
if [ "$with_cooling" -eq 1 ]; then
   echo "Plotting phase space diagram"
   python3 rhoTPlot.py -i 10 -s 'snap/snapshot'
   python3 rhoTPlot.py -i 20 -s 'snap/snapshot'
   python3 rhoTPlot.py -i 30 -s 'snap/snapshot'
   python3 rhoTPlot.py -i 40 -s 'snap/snapshot'
   python3 rhoTPlot.py -i 50 -s 'snap/snapshot'
   python3 rhoTPlot.py -i 100 -s 'snap/snapshot'

   # Movie time!
   python3 rhoTPlot.py -i 0 -f 100 -s 'snap/snapshot'
fi

echo "Plotting radial momentum profile"
python3 radial_momentum_profile.py snap/snapshot_0010.hdf5 \
	-o radial_momentum_profile_0010.png
python3 radial_momentum_profile.py snap/snapshot_0020.hdf5 \
	-o radial_momentum_profile_0020.png
python3 radial_momentum_profile.py snap/snapshot_0030.hdf5 \
	-o radial_momentum_profile_0030.png
python3 radial_momentum_profile.py snap/snapshot_0040.hdf5 \
	-o radial_momentum_profile_0040.png
python3 radial_momentum_profile.py snap/snapshot_0050.hdf5 \
	-o radial_momentum_profile_0050.png
python3 radial_momentum_profile.py snap/snapshot_0100.hdf5 \
	-o radial_momentum_profile_0100.png

# Plot the Sedov solution
echo "Plotting Sedov solution fit"
for ((i=10; i<=100; i+=10)); do
  SNAP_NUMBER=$(printf "%04d" $i)
  INPUT_FILE="snap/snapshot_${SNAP_NUMBER}.hdf5"
  OUTPUT_FILE="sedov_solution_${SNAP_NUMBER}.png"

  echo "Processing ${INPUT_FILE} and saving to ${OUTPUT_FILE}..."

  python3 plot_sedov_solution.py "${INPUT_FILE}" -o "${OUTPUT_FILE}" --step 0.01

  echo "Done with snapshot ${i}."
done

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
	mv *.png $run_name
	mv *.mp4 $run_name
    fi
fi
