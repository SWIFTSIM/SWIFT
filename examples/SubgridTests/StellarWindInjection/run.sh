#!/bin/bash

################# Variation parameters #################
n_threads=${n_threads:=8}  # Number of threads to use
star_type=( 0 ) # Number corrresponding to stellar type : 0=single, 1=continuous_IMF, 2=entire_IMF
star_type_name=( single ) # The name the repository will take for the stellar type
coeff=( 1 ) # The coefficient of pre-supernovae feedback
star_mass=${star_mass:=16} # In solar mass
Z=${Z:=2e-3} # Metallicity (not in solar metallicity)

################## makeIC constant parameters ##################
level=${level:=6}  # Number of particles = 2^(3*level)
jeans_length=${jeans_length:=1}  # Jeans wavelenght in unit of the boxsize
gas_density=${gas_density:=0.1} # Gas density in atom/cm^3
gas_particle_mass=${gas_particle_mass:=10} # Mass of the gas particles
num_star=${num_star:=1} # number of stars

set -e 

# Get the Grackle cooling table
if [ ! -e CloudyData_UVB=HM2012_high_density.h5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    ./getGrackleCoolingTable.sh
fi

if [ ! -e POPII.hdf5 ] || [ ! -e POPIII_PISNe ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    ./getChemistryTable.sh
fi

echo "The star masse = ${star_mass[@]}"
echo "The metallicitie = ${Z[@]}"
echo "star_type = $star_type"
echo "coeff = $coeff"

##################### main script #####################

# Define the directory name based on the parameters
str="test_Star_Type_${star_type_name}_Star_Mass_${star_mass}_Z_${Z}_Coeff_${coeff}"

if [[ "$1" == "SN" ]]; then
    str="${str}_with_SN"
fi

# Create directory if not already present
if [ -d "$str" ]; then
    echo "$str directory exists. Its content will be removed."
    rm -r $str
else
    echo "$str directory does not exist. It will be created."
fi
mkdir $str
cd $str

# Continuous stars are multiplied by a factor of 100, as well as gas particle
mutliply=$(echo "$star_type 0" | awk '{ print ($1 != $2) ? "true" : "false" }')
if [[ "$mutliply" == "true" ]]; then
    star_mass=$( echo "$star_mass 100" | awk '{print $1 * $2 }' )
    gas_particle_mass=$( echo "$gas_particle_mass 100" | awk '{print $1 * $2 }' )
fi

# Generate the Initial Conditions
echo "Generating initial conditions to run the example..."
python3 ../makeIC.py --level $level --lJ $jeans_length --rho $gas_density \
    --mass $gas_particle_mass --star_mass $star_mass \
    --star_type $star_type --num_star $num_star \
    -o ICs_homogeneous_box.hdf5

# Start the simulation
echo "Starting simulation"
if [[ "$2" == "SN" ]]; then
    ../../../../swift --hydro --sync --limiter --external-gravity --stars --sinks --feedback --threads=$n_threads --param="GEARChemistry:initial_metallicity:$Z" --param="GEARFeedback:supernovae_efficiency:0.1" --param="GEARFeedback:pre_supernovae_efficiency:$coeff" ../params.yml 2>&1 | tee output.log
else 
    ../../../../swift --hydro --sync --limiter --external-gravity --stars --sinks --feedback --threads=$n_threads --param="GEARChemistry:initial_metallicity:$Z" --param="GEARFeedback:pre_supernovae_efficiency:$coeff" ../params.yml 2>&1 | tee output.log
fi
cd ..
python3 verification.py --verbose "$str/output.log" | tee verification.txt


