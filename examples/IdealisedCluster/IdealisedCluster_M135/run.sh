#!/bin/bash

with_subgrid={$with_subgrid="EAGLE"} # EAGLE or GEAR

if [ ! -e H135_Tmin65.hdf5 ]
then
    echo "Fetching initial conditions for the idealised cluster example..."
    ./getIC.sh
fi

if [ "$with_subgrid" = "EAGLE" ]
then
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
else

    scripts_location="../../GEAR_ICs_and_SCRIPTS"

    # Get the Grackle cooling table
    if [ ! -e CloudyData_UVB=HM2012.h5 ]; then
	echo "Fetching the Cloudy tables required by Grackle..."
	$scripts_location/getGrackleCoolingTable.sh
    fi

    if [ ! -e POPIIsw.h5 ]; then
	echo "Fetching the chemistry tables..."
	$scripts_location/getChemistryTable.sh
    fi
fi

../../../swift --threads=16 --feedback --external-gravity --self-gravity --stars --star-formation --cooling --temperature --hydro --limiter --sync --black-holes idealised_cluster_M135.yml 2>&1 | tee output.log
