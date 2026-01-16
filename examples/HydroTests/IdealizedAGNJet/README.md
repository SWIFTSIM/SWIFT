An idealized AGN jet hydro test. A conical reservoir of particles is placed in the centre, and particles are kicked radially outwards depending on a jet power.

The particles are treated as fixed boundary particles, until enough time has elapsed that, given a particle's ID, it should be launched. It is kicked, and its ID is increased by 1e7 (arbitary large number), so it is no longer treated as a fixed boundary particle.

The code should be configured as follows to run the example: ./configure --with-hydro=sphenix --with-forcing=idealized-agn-jet

To run, simply execute or submit run.sh, which will create the ICs using the script write_ICs.py, and also create a few analysis plots after the simulation has finished. These will end up in the analysis_scripts directory.

When changing any parameters, please make sure they are changed in both the ICs script and the parameter file, if necessary.
