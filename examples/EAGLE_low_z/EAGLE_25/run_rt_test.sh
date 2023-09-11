#!/bin/bash

# Run a RT test with EAGLE ICs.
#
# To use GEAR-RT, configure SWIFT with
#
# --with-stars=basic --with-hydro=gizmo-mfv --with-riemann-solver=hllc --with-rt=GEAR_1 --with-rt-riemann-solver=GLF --with-feedback=none
# [technically, any other feedback scheme should work as well.]
#
#
# To use the DEBUG RT scheme, configure SWIFT with
# --with-stars=basic --with-rt=debug --with-feedback=none
# [technically, any other feedback scheme should work as well.]


# Generate the initial conditions if they are not present.
if [ ! -e EAGLE_ICs_25.hdf5 ]
then
    echo "Fetching initial conditions for the EAGLE 12Mpc example..."
    ./getIC.sh
fi

../../../swift \
    --hydro --threads=16 --stars --self-gravity \
    --feedback --radiation \
    eagle_25_rt_test.yml

