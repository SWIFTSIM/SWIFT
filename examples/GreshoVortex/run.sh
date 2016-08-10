#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e glass_128.hdf5 ]
then
    echo "Fetching initial glass file for the Gresho-Chan vortex example..."
    ./getGlass.sh
fi
if [ ! -e greshoVortex.hdf5 ]
then
    echo "Generating initial conditions for the Gresho-Chan vortex example..."
    python makeIC.py
fi

# Run SWIFT
../swift -s -t 1 gresho.yml

# Plot the solution
for i in {0..100}
do
    python plotSolution.py $i
done

# Make a movie
mencoder mf://*.png -mf w=645:h=645:fps=12:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o gresho.avi
