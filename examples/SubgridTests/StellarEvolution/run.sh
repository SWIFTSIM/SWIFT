#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e glassCube_32.hdf5 ]
then
    echo "Fetching initial glass file for the Supernovae feedback example..."
    ./getGlass.sh
fi
if [ ! -e stellar_evolution.hdf5 ]
then
    echo "Generating initial conditions for the 3D stellar evolution example..."
    python makeIC.py
fi

# Get the Yield tables
if [ ! -e yieldtables ]
then
    echo "Fetching Yield tables..."
    ./getEagleYieldTable.sh
fi

# Get the solutions
if [ ! -e StellarEvolutionSolution ]
then
    echo "Fetching solutions ..."
    ./getSolutions.sh
fi

../../swift  --feedback --stars --hydro --external-gravity --threads=4 stellar_evolution.yml -P EAGLEChemistry:init_abundance_metal:0.08 2>&1 | tee output.log

python plot_box_evolution.py

../../swift  --feedback --stars --hydro --external-gravity --threads=4 stellar_evolution.yml -P EAGLEChemistry:init_abundance_metal:0.04 2>&1 | tee output.log

python plot_box_evolution.py

../../swift  --feedback --stars --hydro --external-gravity --threads=4 stellar_evolution.yml -P EAGLEChemistry:init_abundance_metal:0.01 2>&1 | tee output.log

python plot_box_evolution.py

../../swift  --feedback --stars --hydro --external-gravity --threads=4 stellar_evolution.yml -P EAGLEChemistry:init_abundance_metal:0.001 2>&1 | tee output.log

python plot_box_evolution.py

../../swift  --feedback --stars --hydro --external-gravity --threads=4 stellar_evolution.yml -P EAGLEChemistry:init_abundance_metal:0.0001 2>&1 | tee output.log

python plot_box_evolution.py
