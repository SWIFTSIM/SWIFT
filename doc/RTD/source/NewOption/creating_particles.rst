.. Adding new schemes
   Darwin Roduit, 5th January 2025

.. _new_option_creating_new_particles:

On-the-fly particle creation
----------------------------

Creating particles on the fly introduces a set of subtleties that need to be addressed. You can find a working example in ``star_formation`` and ``star_formation_sink``, where gas or sink particles spawn stars. *Particle creation* means that we are changing the total number of particles during the simulation.

Particle creation differs from particle type conversion, where the total number of particles remains the same, but the number per type can change. You can find an example of particle conversion within the star formation module.

Updating the maximal displacement
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Swift does not rebuild the cell hierarchy at every step, as it would be computationally inefficient. Instead, we rebuild it when it is not sufficiently good. However, between two tree rebuilds, the particles can move out of their leaf cells, i.e. the cell at the bottom of the tree hierarchy. It's not a problem provided that we know by how much particles have moved between two tree rebuilds, which is called particle *displacement* or *offset*. The displacement information is stored in ``X.dx_max_part``, where ``X`` is the particle type. If this particle type has an associating sorting operation, the particle also owns a ``X.dx_max_sort`` that may be updated if no sort operation are scheduled after spawning new particles. The cells also store the maximum displacement value among the particles they contain.

When we create a new particle, we need to **update the maximal displacement** with respect to its leaf cell. The function ``cell_spawn_new_X_from_Y()`` updates these variables using the parent's particle displacement, as well as setting the newborn particle's position to the parent's one. However, when spawning a new particle, it might be desirable to give a different position than the parent's. In that case, we need to update ``X.dx_max_part`` and ``X.dx_max_sort``. The cells' value must be updated as well.

For star creation, the displacement update per particle is handled  in ``runner_do_star_formation_sink()`` and ``runner_do_star_formation()`` in ``runner_others.c``. Note that ``spart.dx_max_sort`` is not updated as a sort is performed after forming stars. For the cells, a similar machinery is implemented after creating a new particle. Then, the cell's displacement needs to be propagated to the cell hierarchy.
