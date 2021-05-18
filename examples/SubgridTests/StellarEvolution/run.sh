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

if [ ! -e photometry ]
then
    echo "Fetching EAGLE photometry tables..."
    ../getEaglePhotometryTable.sh
fi

# Get the solutions
if [ ! -e StellarEvolutionSolution ]
then
    echo "Fetching solutions ..."
    ./getSolutions.sh
fi

../../swift  --temperature --feedback --stars --hydro --external-gravity --threads=4 stellar_evolution.yml -P EAGLEChemistry:init_abundance_metal:0.08 -P EAGLEChemistry:init_abundance_Hydrogen:0.71 -P EAGLEChemistry:init_abundance_Helium:0.21 2>&1 | tee output_0p08.log

python3 plot_box_evolution.py

../../swift  --temperature --feedback --stars --hydro --external-gravity --threads=4 stellar_evolution.yml -P EAGLEChemistry:init_abundance_metal:0.04 -P EAGLEChemistry:init_abundance_Hydrogen:0.74 -P EAGLEChemistry:init_abundance_Helium:0.23 2>&1 | tee output_0p04.log

python3 plot_box_evolution.py

../../swift  --temperature --feedback --stars --hydro --external-gravity --threads=4 stellar_evolution.yml -P EAGLEChemistry:init_abundance_metal:0.01 2>&1 | tee output_0p01.log

python3 plot_box_evolution.py

../../swift  --temperature --feedback --stars --hydro --external-gravity --threads=4 stellar_evolution.yml -P EAGLEChemistry:init_abundance_metal:0.001 2>&1 | tee output_0p001.log

python3 plot_box_evolution.py

../../swift  --temperature --feedback --stars --hydro --external-gravity --threads=4 stellar_evolution.yml -P EAGLEChemistry:init_abundance_metal:0.0001 2>&1 | tee output_0p0001.log

python3 plot_box_evolution.py
