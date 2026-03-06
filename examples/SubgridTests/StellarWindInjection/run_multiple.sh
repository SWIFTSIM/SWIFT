#!/bin/bash
################# Variation parameters #################
n_threads=${n_threads:=12} # Number of threads to use
star_type=( 0 ) # Number corrresponding to stellar type : 0=single, 1=continuous_IMF, 2=entire_IMF
star_type_name=( single continuous ) # The name the repository will take for the stellar type
coeff=( 1 ) # The coefficient of pre-supernovae feedback

# Logarithmic spacing
star_mass_min=${star_mass_min:=0.5}
star_mass_max=${star_mass_max:=400}
star_mass_count=${star_mass_count:=2}
Z_min=${Z_min:=1e-11}
Z_max=${Z_max:=2e-1}
Z_count=${Z_count:=2}

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

if [ ! -e POPII.hdf5 ] || [! -e POPIII_PISNe.hdf5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    ./getChemistryTable.sh
fi

star_mass=($(awk -v min=$star_mass_min -v max=$star_mass_max -v count=$star_mass_count 'BEGIN{for(i=0;i<count;i++)print min*(max/min)^(i/(count-1))}'))
Z=($(awk -v min=$Z_min -v max=$Z_max -v count=$Z_count 'BEGIN{for(i=0;i<count;i++)print min*(max/min)^(i/(count-1))}'))

echo "The star masses are : ${star_mass[@]}"
echo "The metallicities are : ${Z[@]}"

##################### main script #####################
# Calculate the total number of combinations
echo "Debug : ${#Z[@]} ${#star_type[@]} ${#star_mass[@]} ${#coeff[@]}"
total_combinations=$((${#Z[@]} * ${#star_type[@]} * ${#star_mass[@]} * ${#coeff[@]}))
echo "total combinations : ${total_combinations}"

# Determine the combination index based on SLURM_ARRAY_TASK_ID
for ((index=0; index<total_combinations; index++)); do 

    echo "SLURM_ARRAY_TASK_ID: $index"

    # Calculate the correct set of parameters for the current job
    z_index=$((index % ${#Z[@]}))
    star_type_index=$(( (index / ${#Z[@]}) % ${#star_type[@]} ))
    star_mass_index=$(( (index / (${#Z[@]} * ${#star_type[@]})) % ${#star_mass[@]} ))
    coeff_index=$(( (index / (${#Z[@]} * ${#star_type[@]} * ${#star_mass[@]})) % ${#coeff[@]} ))

    echo "z_index: $z_index"
    echo "star_type_index: $star_type_index"
    echo "star_mass_index: $star_mass_index"
    echo "coeff_index: $coeff_index"

    # Get the values for the parameters
    z=${Z[$z_index]}
    star_type_value=${star_type[$star_type_index]}
    star_type_name_value=${star_type_name[$star_type_index]}
    star_mass_value=${star_mass[$star_mass_index]}
    coeff_value=${coeff[$coeff_index]}

    # Define the directory name based on the parameters
    str="test_Star_Type_${star_type_name_value}_Star_Mass_${star_mass_value}_Z_${z}_Coeff_${coeff_value}"

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

    mutliply=$(echo "$star_type_value 0" | awk '{ print ($1 != $2) ? "true" : "false" }')
    #echo "$mutliply"
    if [[ "$mutliply" == "true" ]]; then
        star_mass_value=$( echo "$star_mass_value 100" | awk '{print $1 * $2 }' )
        gas_particle_mass=$( echo "$gas_particle_mass 100" | awk '{print $1 * $2 }' )
        #echo "should be * 100"
    fi

    echo "Star type: $star_type_value -> $star_type_name_value"
    echo "Star mass: $star_mass_value"
    echo "Metallicity: $z"
    echo "Stellar wind efficiency: $coeff_value"

    # Generate the Initial Conditions
    echo "Generating initial conditions to run the example..."
    python3 ../makeIC.py --level $level --lJ $jeans_length --rho $gas_density \
        --mass $gas_particle_mass --star_mass $star_mass_value \
        --star_type $star_type_value --num_star $num_star \
        -o ICs_homogeneous_box.hdf5

    # Start the simulation
    echo "Starting simulation"
    if [[ "$2" == "SN" ]]; then
        ../../../../swift --hydro --sync --limiter --external-gravity --stars --sinks --feedback --threads=$n_threads --param="GEARChemistry:initial_metallicity:$z" --param="GEARFeedback:supernovae_efficiency:0.1" --param="GEARFeedback:pre_supernovae_efficiency:$coeff_value" ../params.yml 2>&1 | tee output.log
    else 
        ../../../../swift --hydro --sync --limiter --external-gravity --stars --sinks --feedback --threads=$n_threads --param="GEARChemistry:initial_metallicity:$z" --param="GEARFeedback:pre_supernovae_efficiency:$coeff_value" ../params.yml 2>&1 | tee output.log
    fi
    cd ..
done

python3 ./verification.py --verbose ./test_*/output.log | tee verification.txt
