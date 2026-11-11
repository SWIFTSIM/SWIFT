Orszag-Tang Vortex MHD test in 3D.
This test requires setting mu_0 = 1.
We recommend running with a quartic spline kernel ;
the parameters in the .yml file are chosen to work best	in that	configuration.
The fiducial set up uses stacked BCC lattices in the ICs, 
but you can instead use glass files by running 

./getGlass.sh

and addapting the makeIC.py script accordingly.
