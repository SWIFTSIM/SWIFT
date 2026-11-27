Blob Test
=========

The test in this folder is identical to the one presented in Braspenning+ 2022
(https://ui.adsabs.harvard.edu/abs/2022arXiv220313915B/abstract)
and can be used to reproduce those results.

You can run this test with any of the hydro schemes supported by SWIFT; no
special configure option is needed.

Initial conditions
------------------

Download the ICs with `./getICs.sh`. It fetches the eight files
`blob_<N>_<chi>.hdf5` (only those not already present), where

+ `<N>` (16, 32, 64 or 128) is the resolution: the wind has N^3 particles per
  unit cube, so the 4:1:1 box contains about 4 N^3 particles;
+ `<chi>` (10 or 100) is the density contrast between the cloud and the wind.

`blob.yml` reads `blob_16_100.hdf5`. To use another set, change
`InitialConditions:file_name` in `blob.yml` or pass it on the command line.

Running
-------

`plot_cloud_evolution.py` expects the snapshots of a run in a sub-directory
named after the resolution (e.g. `./128/blob_0000.hdf5`), so run each
resolution in its own directory. For instance, from this folder:

    ./getICs.sh
    mkdir -p 128 && cd 128
    ../../../../swift --pin --hydro --limiter --threads=28 \
        --param=InitialConditions:file_name:../blob_128_10.hdf5 ../blob.yml
    cd ..

Plotting
--------

`python3 plot_cloud_evolution.py` makes Fig. 7 of Braspenning+ 2022. The
resolution and density contrast to analyse are set at the top of the script
(`resolution = 128` and `chi = 10` by default) and must match the run. The
figure (`cloud_evolution.png`) shows:

+ top row: the evolution of the mass of dense gas and of
  intermediate-temperature gas;
+ second and third rows: the same quantities after resampling the gas onto
  16^3 and 32^3 grids.
