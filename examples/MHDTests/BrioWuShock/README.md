Brio & Wu 3D MHD shock tube.
This test requires setting mu_0 = 1 and gamma = 2.
We recommend running with a quintic spline kernel ;
the parameters in the .yml file are chosen to work best	in that	configuration.
The fiducial set up uses stacked BCC lattices in the ICs, 
but you can instead use glass files by running 

./getGlass.sh

and addapting the makeIC.py script accordingly.