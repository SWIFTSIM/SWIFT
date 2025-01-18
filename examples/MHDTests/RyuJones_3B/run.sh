
#!/bin/bash

# Generate the initial conditions file if not present.
if [ ! -f ./BCCglassCube_32.hdf5 ]
then
    echo "Generating a BCC unit cell, copies of which are to be stacked to generate the Ryu & Jones 3B shock tube ..."
    python3 makeBCC.py -n 32
fi
if [ ! -f ./RyuJones_3B.hdf5 ]
then
    echo "Generating initial conditions for the Ryu & Jones 3B shock tube example..."
    python3 makeIC.py -n 4
fi

# Run the example with SWIFT
../../../swift --hydro --threads=4 RyuJones_3B.yml 2>&1 | tee output.log

# Plot the calculated solution at time t=0.1
python3 plotSolution.py RyuJones_3B_0011.hdf5 RyuJones_3B_0011.png
