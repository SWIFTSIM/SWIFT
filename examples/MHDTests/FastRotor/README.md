Fast Rotor MHD test in 3D, following :
Lodrillo P. & Del Zanna L., 2000
Stasyszyn et al., 2013
Seo J. & Ryu D., 2023

This test requires setting mu_0 = 1 and gamma = 7/5.
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

