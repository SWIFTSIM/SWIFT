#!/bin/bash

printf "Running simulation..."
python3 redo_ICs_CDM_PID.py # fix PID=0 issue
../../../swift --self-gravity --threads=14 params_cdm.yml 2>&1 | tee output.log