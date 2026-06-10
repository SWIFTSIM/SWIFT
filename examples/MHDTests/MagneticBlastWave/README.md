Magnetic Blast Wave test in 3D, following :
Pakmor R., Bauer A. & Springel V., 2011 
Hopking P. F. & Raives M. J., 2016

This test requires setting mu_0 = 1.
We recommend running with a quartic spline kernel ;
the parameters in the .yml file are chosen to work best	in that	configuration.

It is also possible to automatically compile and run the test for all available
MHD schemes, sequentially (files in the directoy with the _schemes suffix serve
this purpose). To do so, run :

./run_schemes.sh [what Scheme] [FOLDER_TAIL]

where :

[what scheme]: 
vep: vector potentials
odi: Oresti's direct induction
fdi: simple direct induction
[FOLDER_TAIL]: 
the trailing name for the folders created
