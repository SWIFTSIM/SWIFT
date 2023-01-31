.. Light Cones
   John Helly 29th April 2021

.. _lightcone_algorithm_description_label:

Light Cone Output Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In cosmological simulations it is possible to specify the location of
an observer in the simulation box and have SWIFT output information
about particles in the simulation as they cross the observer's past
light cone.

Whenever a particle is drifted the code checks if any periodic copy of
the particle crosses the lightcone during the drift, and if so that
copy of the particle is buffered for output. As an optimization, at the
start of each time step the code computes which periodic copies of the
simulation box could contribute to the light cone and only those copies
are searched. When drifting the particles in a particular cell the list of
replications is further narrowed down using the spatial extent of the
cell.

Particles can be output directly to HDF5 files or accumulated to healpix
maps corresponding to spherical shells centred on the observer.



