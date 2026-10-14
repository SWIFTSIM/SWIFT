Blob Test
=========

The test in this folder is identical to the one presented in Braspenning+ 2022
(https://ui.adsabs.harvard.edu/abs/2022arXiv220313915B/abstract)
and can be used to reproduce those results

You can use the provided ICs, where the first value indicates the particle number
and the second the density contrast. Switching ICs requires indicating the new ICs
in the parameter file 'blob.yml'
You can run this test with any of the hydro schemes supported by SWIFT.

Figures similar to those in Braspenning+ 2022 can be made using 'plot_cloud_evolution.py',
this produces:
+ Evolution of the mass of dense gas
+ Evolution of the mass of intermediate-temperature gas
The second and third row show the same evolutions by downsampled to a lower resolution grid.
This is identical to Fig. 7 in Braspenning+ 2022

The simulation can be run using this command:
./swift --pin --hydro --limiter --threads=28 blob.yml
