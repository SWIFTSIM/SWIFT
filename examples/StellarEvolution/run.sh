#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e stellar_evolution.hdf5 ]
then
    echo "Generating initial conditions for the 3D stellar evolution example..."
    python makeIC.py
fi

# Run mass enrichment check
#../swift --limiter --feedback --stars --hydro --external-gravity --threads=4 --param=Stars:energy_testing:0 --param=TimeIntegration:time_end:1e-4 --param=Stars:feedback_timescale:1e-4 stellar_evolution.yml 2>&1 | tee output_enrichment.log
#
#python check_stellar_evolution.py
#
## Run continuous heating check
#../swift --limiter --feedback --stars --hydro --external-gravity --threads=4 --param=Stars:energy_testing:1 --param=Stars:continuous_heating:1 stellar_evolution.yml 2>&1 | tee output_continuous.log
#
#python check_continuous_heating.py
#
## Run stochastic check
#../swift --limiter --feedback --stars --hydro --external-gravity --threads=4 --param=Stars:energy_testing:1 stellar_evolution.yml 2>&1 | tee output_stochastic_1.log
#
#python check_stochastic_heating.py

../swift --limiter --feedback --stars --hydro --external-gravity --threads=4 --param=Stars:energy_testing:1 stellar_evolution.yml 2>&1 | tee output.log

python check_stochastic_heating.py
#
#../swift --limiter --feedback --stars --hydro --external-gravity --threads=4 --param=Stars:energy_testing:1  stellar_evolution.yml 2>&1 | tee output_stochastic_3.log
#
#python check_stochastic_heating.py
#
#../swift --limiter --feedback --stars --hydro --external-gravity --threads=4 --param=Stars:energy_testing:1  stellar_evolution.yml 2>&1 | tee output_stochastic_4.log
#
#python check_stochastic_heating.py
