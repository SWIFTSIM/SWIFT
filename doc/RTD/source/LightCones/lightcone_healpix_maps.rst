.. Light Cones
   John Helly 29th April 2021

.. _lightcone_healpix_maps_label:

Light Cone HEALPix Maps
~~~~~~~~~~~~~~~~~~~~~~~

SWIFT can accumulate particle properties to HEALPix maps as they
cross the observer's past light cone. Each map corresponds to a
spherical shell centred on the observer. When a particle crosses
the lightcone its distance from the observer is calculated and the
particle's contribution is added to a buffer so that at the end of
the time step it can be added to the corresponding HEALPix map.

Maps can be generated for multiple concentric shells and multiple
quantities can be accumulated for each shell. The HEALPix map for a
shell is allocated and zeroed out when the simulation first reaches
a redshift where particles could contribute to that map. The map is
written out and deallocated when the simulation advances to a point
where there can be no further contributions. In MPI runs the pixel
data for the maps are distributed across all MPI ranks.

Updates to the maps are buffered in order to avoid the need for
communication during the time step. At the end of the step if any
MPI rank has a large amount of updates buffered then all pending
updates will be applied to the pixel data.

For gas particles, the HEALPix maps are smoothed using a projected
version of the same kernel used for the hydro calculations. Other
particle types are not smoothed.

The code writes one output file for each spherical shell. In MPI mode
all ranks write to the same file using parallel HDF5. If maps of
multiple quantities are being made they will be written to a single
file as separate 1D datasets with one element per pixel.
