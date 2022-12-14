.. Orphan Particle Output
   John Helly 18th September 2020

Orphan Particle Output
======================

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Contents: 

When running VELOCIraptor on the fly SWIFT can write additional output
files which contain the positions, velocities and IDs of all particles
which have ever been the most bound particle in a halo. The next section
describes the motivation for this feature. Subsequent sections describe how
to enable it and the output format.

Orphan galaxies in semi-analytic galaxy formation models
--------------------------------------------------------

SWIFT and VELOCIraptor output may be used as the basis for semi-analytic
galaxy formation models. A SWIFT simulation of a cosmological box is run
with on-the-fly VELOCIraptor output, halo merger trees are constructed
and a semi-analytic modelling code is used to populate the halos with
galaxies. These models often assume that satellite galaxies in a halo 
continue to survive after their associated dark matter subhalo is no
longer identified by the halo finder. In that case the galaxy (referred to
as an 'orphan' because it has no parent subhalo) may be associated with
the most bound particle of the subhalo at the last time it was identified.

Enabling orphan particle output
-------------------------------

Orphan particle output can only be enabled when running VELOCIraptor on
the fly. It is enabled by running configure with the flag
``--enable-velociraptor-orphans``. E.g.

.. code-block:: bash

  ./configure --with-velociraptor=/path/to/VELOCIraptor-STF/src --enable-velociraptor-orphans

If SWIFT is configured with this flag then whenever VELOCIraptor is invoked
during a simulation additional output files will be generated alongside the
usual VELOCIraptor output.

Output format
-------------

The output files are generated in the same location and with the same base
name as the VELOCIraptor output but with a .orphans suffix. If parallel
HDF5 is available then there will be one output file per VELOCIraptor
invocation. Otherwise there is one output file for each MPI task. Non-MPI
runs always generate one file per VELOCIraptor invocation.

The output files are written as HDF5 and are in a similar format to the
:ref:`snapshots`, including the meta-data groups. Note that all particles
are written to the PartType1 group irrespective of their type. This feature
is primarily intended for use with dark matter only cosmological simulations.

The PartType1 group always contains the datasets Coordinates, ParticleIDs
and Velocities in the same units as in the snapshot files and with the same
meta-data attributes attached.
