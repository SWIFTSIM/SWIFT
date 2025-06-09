#!/bin/bash

if [ ! -e fid.hdf5 ] 
then     
    echo "Fetching initial conditions for the isolated galaxy example..."
    ./getIC.sh
fi

if [ ! -e coolingtables ] 
then     
    echo "Fetching EAGLE cooling tables for the isolated galaxy example..."
    ./getEagleCoolingTable.sh
fi

if [ ! -e UV_dust1_CR1_G1_shield1.hdf5 ]
then
    echo "Fetching PS2020 cooling tables for the isolated galaxy example..."
    ./getPS2020CoolingTables.sh
fi

if [ ! -e yieldtables ] 
then     
    echo "Fetching EAGLE stellar yield tables for the isolated galaxy example..."
    ./getYieldTable.sh
fi

if [ ! -e photometry ]
then
    echo "Fetching EAGLE photometry tables..."
    ../getEaglePhotometryTable.sh
fi

# Adding magnetic fields.
echo "Adding magnetic fields..."
python3 makeIC.py fid.hdf5 fid_B.hdf5 

# Run SWIFT with external potential, cooling, self-gravity and star formation
#../../../swift --external-gravity --self-gravity --hydro --cooling --stars --star-formation --threads=4 isolated_galaxy.yml 2>&1 | tee output.log

# Run SWIFT with external potential, cooling, self-gravity, star formation and feedback
../../../swift --external-gravity --self-gravity --hydro --cooling --stars --star-formation --feedback --sync --threads=4 isolated_galaxy.yml 2>&1 | tee output.log
