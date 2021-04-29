.. Light Cones
   John Helly 29th April 2021

.. _lightcone_running_label:

Running SWIFT with Light Cone Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To produce light cone particle output swift must be configured
with ``--enable-lightcone``. Additionally, making HEALPix maps
requires the HEALPix C library. If using MPI then parallel HDF5
is also required.

Several parameters must be added to the SWIFT parameter file. For
example:

```
Lightcone:

  # Common parameters
  basename: lightcone
  observer_position: [35.5, 78.12, 12.45]
  buffer_chunk_size: 10000

  # Particle output parameters
  z_min_for_particles:    0.0
  z_max_for_particles:    0.05
  use_gas:           1
  use_dm:            1
  use_dm_background: 0
  use_stars:         0
  use_black_hole:    0
  use_neutrino:      0
  max_particles_buffered: 10000
  hdf5_chunk_size:        10000

  # Healpix map parameters
  nside:                512
  radius_file:          ../../shell_radii.txt
  max_updates_buffered: 100000
  map_names:            [TotalMass]
```
