#!/bin/bash

################# Variation parameters #################
n_threads=${n_threads:=8}  # Number of threads to use
star_type=( 0 ) # Number corrresponding to stellar type : 0=single, 1=continuous_IMF, 2=entire_IMF
star_type_name=( single ) # The name the repository will take for the stellar type
coeff=( 1 ) # The coefficient of pre-supernovae feedback
star_mass=${star_mass:=32} # In solar mass
Z=${Z:=2e-3} # Metallicity (not in solar metallicity)
gas_particle_mass=${gas_particle_mass:=76} # Mass of the gas particles

################## makeIC constant parameters ##################
level=${level:=6}  # Number of particles = 2^(3*level)
jeans_length=${jeans_length:=1}  # Jeans wavelenght in unit of the boxsize
gas_density=${gas_density:=100} # Gas density in atom/cm^3
num_star=${num_star:=1} # number of stars

set -e 

# Get the Grackle cooling table
if [ ! -e CloudyData_UVB=HM2012_high_density.h5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    ./getGrackleCoolingTable.sh
fi

if [ ! -e POPII.hdf5 ] || [ ! -e POPIII_PISNe.hdf5 ]
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
str="${star_type_name}_level${level}_M${star_mass}_Z${Z}_Rho${gas_density}_GasM${gas_particle_mass}"

# Create directory if not already present
if [ -d "$str" ]; then
    echo "$str directory exists. Its content will be removed."
    rm -r $str
else
    echo "$str directory does not exist. It will be created."
fi
mkdir $str
cd $str

# Generate the Initial Conditions
echo "Generating initial conditions to run the example..."
python3 ../makeIC.py --level $level --lJ $jeans_length --rho $gas_density \
    --mass $gas_particle_mass --star_mass $star_mass \
    --star_type $star_type --num_star $num_star \
    -o ICs_homogeneous_box.hdf5

#--------------------- Run the simulation | no cooling part --------------------#
mkdir noCooling
cd noCooling
# Start the simulation
../../../../../swift --hydro --sync --limiter --external-gravity --stars --sinks --feedback  --threads=$n_threads --param="GEARChemistry:initial_metallicity:$Z" --param="GEARFeedback:pre_supernovae_efficiency:$coeff" ../../params.yml 2>&1 | tee output.log

# Generate the plots
python3 ../../gas_profile_analysis.py "./" | tee gas_profile_analysis.txt

python3 ../../bubble_evolution.py --log "radial_peak_positions.txt"
cd ..
#--------------------- Run the simulation | cooling part --------------------#
mkdir withCooling
cd withCooling
# Start the simulation
../../../../../swift --hydro --sync --limiter --external-gravity --cooling --stars --sinks --feedback  --threads=$n_threads --param="GEARChemistry:initial_metallicity:$Z" --param="GEARFeedback:pre_supernovae_efficiency:$coeff" ../../params.yml 2>&1 | tee output.log 
# Generate the plots
python3 ../../gas_profile_analysis.py "./" | tee gas_profile_analysis.txt 
python3 ../../bubble_evolution.py --log "radial_peak_positions.txt"




