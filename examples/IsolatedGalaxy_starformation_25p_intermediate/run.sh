#!/bin/bash

if [ ! -e starforming_intermediate.hdf5 ]
then 
  ./getIC.sh
fi 

../swift --threads=32 --external-gravity --self-gravity --stars --star-formation --cooling --temperature --hydro isolated_galaxy.yml 2>&1 | tee output.log
