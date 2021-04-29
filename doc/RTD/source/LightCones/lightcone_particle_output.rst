.. Light Cones
   John Helly 29th June 2021

.. _lightcone_particle_output_label:

Light Cone Particle Output
~~~~~~~~~~~~~~~~~~~~~~~~~~

SWIFT can output particles to HDF5 output files (similar to the
snapshots) as they cross the observer's light cone. During each time
step, any particles which cross the light cone are added to a buffer.
If this buffer is large at the end of the step then its contents
are written to an output file. In MPI runs each MPI rank writes its
own output file and decides independently when to flush its particle
buffer.

A new output file is started whenever restart files are written. This
allows the code to automatically continue from the point of the restart
dump if the run is interrupted. Any files written after the restart
dump will be overwritten when the simulation is resumed, preventing
duplication of particles in the light cone output.

The output files have names of the form ``basename_XXXX.Y.hdf5``, where
XXXX numbers the files written by a single MPI rank and Y is the index
of the MPI rank.

The output files contain one HDF5 group for each particle type. Within
each group there are datasets corresponding to particle properties in
a similar format to the snapshots.

