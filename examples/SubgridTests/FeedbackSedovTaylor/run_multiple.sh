#!/bin/bash

# Number of different seeds to run
n_runs=${n_runs:=10}  # Change this to the desired number of runs

# Base directory to store results
RESULTS_DIR=${run_name:="simulation_results"}
mkdir -p "$RESULTS_DIR"

# Fixed parameter
export with_cooling=1

# Arrays of gas_particle_mass and level values to iterate over
gas_particle_masses=(10000 1000 100 10)

# Loop over different configurations
for gas_mass in "${gas_particle_masses[@]}"; do
    export gas_particle_mass=$gas_mass
    export boxsize=2 # kpc

    echo "Running simulation with gas_particle_mass=$gas_mass, level=$lvl..."
    
    # Run the simulation script
    ./run.sh
    
    # Move output to a new directory
    RUN_DIR="$RESULTS_DIR/run_mass_${gas_mass}_level_${lvl}"
    mkdir -p "$RUN_DIR"
    mv snap "$RUN_DIR/"
    mv output.log "$RUN_DIR/"
    mv timesteps.txt "$RUN_DIR/"
    mv statistics.txt "$RUN_DIR/"
    mv dependency_graph_0.csv "$RUN_DIR/"
    mv *.png "$RUN_DIR/"
    mv *.mp4 "$RUN_DIR/"
    
    echo "Finished run $i with gas_particle_mass=$gas_mass, level=$lvl. Results stored in $RUN_DIR"
done
