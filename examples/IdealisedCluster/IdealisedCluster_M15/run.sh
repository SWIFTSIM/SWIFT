#!/bin/bash

if [ ! -e M15_fiducial.hdf5 ] 
then     
    echo "Fetching initial conditions for the idealised cluster example..."
    ./getIC.sh
fi

if [ ! -e UV_dust1_CR1_G1_shield1.hdf5 ]
then
    echo "Fetching PS2020 cooling tables for the isolated galaxy example..."
    ../getPS2020CoolingTables.sh
fi

if [ ! -e yieldtables ] 
then     
    echo "Fetching EAGLE stellar yield tables for the isolated galaxy example..."
    ../getYieldTable.sh
fi

if [ ! -e photometry ]
then
    echo "Fetching EAGLE photometry tables..."
    ../getEaglePhotometryTable.sh
fi

../../../swift --threads=16 --feedback --external-gravity --self-gravity --stars --star-formation --cooling --temperature --hydro --limiter --sync --black-holes idealised_cluster_M15.yml 2>&1 | tee output.log
