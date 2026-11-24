#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e primordialBox.hdf5 ]
then
    echo "Generating the initial conditions of the box..."
    python3 makeIC.py
fi

SWIFT=${SWIFT:-../../../swift}

# With --no-cosmic-rays the same two runs go without cosmic rays.
suffix=""
yml3=network3.yml
yml4=network4.yml
if [ "$1" = "--no-cosmic-rays" ]
then
    suffix="_noCR"
    for n in 3 4
    do
        sed -e '/cosmic_ray/d' -e "s/^\(  basename: *network$n\)/\1$suffix/" network$n.yml > network$n$suffix.yml
    done
    yml3=network3$suffix.yml
    yml4=network4$suffix.yml
fi

# Run SWIFT with network 3, then with network 4
$SWIFT --hydro --cooling --threads=4 $yml3
$SWIFT --hydro --cooling --threads=4 $yml4

# Print the difference of the two networks
python3 compareNetworks.py network3$suffix network4$suffix
