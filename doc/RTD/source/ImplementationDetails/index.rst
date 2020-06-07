.. Implementation details
   Loic Hausammann, 2020
   Matthieu Schaller, 2020

.. _implementation_details:

Implementation details
======================

Generating new unique IDs
-------------------------

When spawning new particles (not converting them) for star formation or other
similar processes, the code needs to create new unique particle IDs. This is
implemented in the file ``space_unique_id.c`` and can be switched on/off in the
star formation file ``star_formation_struct.h`` by setting the variable
``star_formation_need_unique_id`` to 1 or 0.

The generation of new IDs is done by computing the maximal ID present in the
initial condition (across all particle types) and then attributing two batches
of new, unused IDs to each MPI rank.  The size of each batch is computed in the
same way as the count of extra particles in order to ensure that we will have
enough available IDs between two tree rebuilds (where the extra particles are
regenerated).

When a new ID is requested, the next available ID in the first batch is
returned. If the last available ID in the first batch is requested, we switch to
the next batch of IDs and flag it for regeneration at the next rebuild time.  If
the second batch is also fully used, the code will exit with an error message
[#f1]_. At each tree-rebuild steps, the ranks will request a new batch if
required and make sure the batches are unique across all MPI ranks.

As we are using the maximal ID from the ICs, nothing can be done against the user
providing the maximal integer possible as an ID (that can for instance be the
case in some of the EAGLE ICs as the ID encode their Lagrangian position on a
Peano-Hilbert curve). 


.. [#f1] Thanks to the size of the fresh ID batches, the code should run out of
	 extra particles before reaching this point and triggered a new rebuild
	 if this is allowed by the star formation scheme.
