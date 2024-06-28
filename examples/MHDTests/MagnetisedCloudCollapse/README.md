Magnetised Cloud Collapse

---------------------------------------------------------------------

Test case for investigation of coupling of MHD to gravity.

Generate ICs with makeIC.py. You can provide as a command line
argument the particle number on a side of the cube out of which
we carve the collapsing spherical cloud (all other attributes
of the set up are computed automatically) e.g. :

python3 makeIC.py -n 128 

for a high resolution setup. The script also prints an
indicative value of the gravitational softening length necessary to
resolve the densities of physical interest; addapt accordingly in
the parameter file.

The simulation can then be run as :

../../../swift -s -G --threads=12 magnetised_cloud.yml

Results for snapshot snap.hdf5 can then be plotted and storred in
im.png as :

python3 plotSolution.py snap.hdf5 im.png

The directory also contains a bash script that does all this for you.
