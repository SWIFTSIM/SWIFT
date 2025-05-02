#!/bin/bash

#Creates a directory for the outputs
DIR=output_1 #First test of units conversion
if [ -d "$DIR" ];
then
    echo "$DIR directory exists. Its content will be removed."
    rm $DIR/output_*
else
    echo "$DIR directory does not exists. It will be created."
    mkdir $DIR
fi

DIR=output_2 #Second test of units conversion
if [ -d "$DIR" ];
then
    echo "$DIR directory exists. Its content will be removed."
    rm $DIR/output_*
else
    echo "$DIR directory does not exists. It will be created."
    mkdir $DIR
fi

#Clears the previous figures
echo "Clearing existing figures."
if [ -f "circular_orbits_simulation_kpc.png" ];
then
    rm circular_orbits_simulation_kpc.png
fi

if [ -f "circular_orbits_simulation_Mpc.png" ];
then
    rm circular_orbits_simulation_Mpc.png
fi

if [ -f "deviation_simulation_kpc.png" ];
then
    rm deviation_simulation_kpc.png
fi
if [ -f "deviation_simulation_Mpc.png" ];
then
    rm deviation_simulation_Mpc.png
fi

if [ -f "deviation_from_original_data_simulation_kpc.png" ];
then
    rm deviation_from_original_data_simulation_kpc.png
fi

if [ -f "deviation_from_original_data_simulation_Mpc.png" ];
then
    rm deviation_from_original_data_simulation_Mpc.png
fi

#Clears the IC file
if [ -f "circular_orbits_MW.hdf5" ];
then
    rm circular_orbits_MW.hdf5
fi


#Generates the initial conditions
echo "Generate initial conditions for circular orbits"
if command -v python3 &>/dev/null; then
    python3 makeIC.py
else
    python3 makeIC.py
fi

#Runs the simulation
# self gravity G, external potential g, hydro s, threads t and high verbosity v
echo "Starting the simulation with the first type of units (kpc)... Look at output_1.log for the simulation details."
../../../swift --external-gravity --threads=8 params_unit_1.yml 2>&1 > output_1.log
echo "Simulation ended."

echo "Starting the simulation with the second type of units (Mpc)... Look at output_2.log for the simulation details."
../../../swift --external-gravity --threads=8 params_unit_2.yml 2>&1 > output_2.log
echo "Simulation ended."

#Saves the plots
echo "Save plots of the circular orbits and of the errors"
if command -v python3 &>/dev/null; then
    python3 makePlots.py
else 
    python3 makePlots.py
fi
