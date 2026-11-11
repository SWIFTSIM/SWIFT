Magnetised Cloud Collapse

---------------------------------------------------------------------

Test case for investigation of coupling of MHD to gravity.

The test has been ivestigated with the Direct Induction (Price 2018)
formulation of MHD, and default parameters in this directory have been
chosen to best couple with a Wendland-C2 kernel. Note the need for
a barotropic equation of state and an adiabatic index of 4/3. We
thus suggest you compile SWIFT as :

./configure --with-spmhd=direct-induction --with-kernel=wendland-C2 --with-equation-of-state=barotropic-gas --with-adiabatic-index=4/3

to run this test.

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
