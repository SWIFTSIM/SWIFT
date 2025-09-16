#!/bin/bash

# Number of different seeds to run
n_runs=${n_runs:=10}  # Change this to the desired number of runs
level=${level:=5}  #Number of particles = 2^(3*level)

# Base directory to store results
RESULTS_DIR=${run_name:="simulation_results"}
mkdir -p "$RESULTS_DIR"

export level=$level

# Loop over different seeds
for i in $(seq 1 $n_runs); do
    echo "Running simulation with seed $i..."
    
    # Set the seed for this run
    export seed=$i
    
    # Run the simulation script
    ./run.sh
    
    # Move output to a new directory
    RUN_DIR="$RESULTS_DIR/run_$i"
    mkdir -p "$RUN_DIR"
    mv snap "$RUN_DIR/"
    mv output.log "$RUN_DIR/"
    mv timesteps.txt "$RUN_DIR/"
    mv statistics.txt "$RUN_DIR/"
    mv dependency_graph_0.csv "$RUN_DIR/"
    
    echo "Finished run $i. Results stored in $RUN_DIR"
done
