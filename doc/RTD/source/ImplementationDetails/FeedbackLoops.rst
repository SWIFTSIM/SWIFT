.. Stellar feedback loops
   Darwin Roduit, 2026

.. _stellar_feedback_loops:


Stellar feedback loops
~~~~~~~~~~~~~~~~~~~~~~

In all feedback loops, we only consider stars that are active on the timeline (``spart_is_active()``) and if they are active for feedback (``feedback_is_active()``).  The latter function's behaviour is implementation-dependent and defined for
each feedback module.

Gas particles are considered only if they are not inhibited (``part_is_inhibited()``).

A key architectural detail in SWIFT feedback is that **stars (spart) are always the active particles**. They are the "setters" that initiate the neighbour search. **Gas particles (part) are passive** "getters." Therefore, within the current stellar feedback architecture, the gas particles do not perform any neighbour loops to look for star particles.

The page :ref:`current_dependencies` summarises the different feedback tasks available in SWIFT.

Summary of Feedback Stages
--------------------------

The table below summarises the flow of information across the different tasks:

+---------------+-----------------+-------------------------------------------------------+
| Task Name     | Updated Type    | Purpose                                               |
+===============+=================+=======================================================+
| **Density**   | ``spart``       | Star searches for gas neighbors to determine their    |
|               |                 | smoothing length ``h``.                               |
+---------------+-----------------+-------------------------------------------------------+
| **Prep 1**    | ``part``        | Star updates its neighbors. Used in EAGLE kinetic     |
|               |                 | feedback for XXXXX (TODO).                            |
+---------------+-----------------+-------------------------------------------------------+
| **Prep 2-4**  | ``spart``       | Star updates its own internal data using gas info     |
|               |                 | (e.g., computing normalisation or vector weights).    |
+---------------+-----------------+-------------------------------------------------------+
| **Feedback**  | ``part``        | Star injects physical quantities (energy, momentum,   |
|               |                 | metals) back into the gas neighbours.                 |
+---------------+-----------------+-------------------------------------------------------+

The distinction of which loop updates which particle is relevant, as the **MPI** communication tasks will only exchange the **specified particle type**. Therefore, one cannot update ``spart`` during ``prep1`` since we only communicate the gas ``part``.

Task and Ghost Details
----------------------

While the loops perform the physics, the "ghost" tasks manage the synchronisation
and finalisation of data across cells and MPI ranks.

* **Density Loop**: Computes the basic kernel-weighted properties.
* **Density Ghost**: This task finalises the density computations. It is
  computationally expensive and can benefit from **multiplexing**. (TODO: Add a quick description of multiplexing)
* **Prep 1 Loop**: This is a unique stage where the star particle modifies its gas neighbours *before* the final feedback stage.
* **Prep 2, 3, and 4 Loops**: These are used to perform additional computations. The star reads from the gas to update its own ``feedback_data`` structure.
* **Feedback Loop**: The final stage where the physical interaction occurs. The star uses the accumulated weights from previous loops to distribute, e.g. mass, metals, and momentum to the gas. The  ``part`` and ``xpart`` structures can be modified here.

Implicit Tasks and Ghosts
-------------------------

Several tasks in the feedback graph are "implicit", i.e. they exist in the scheduler but perform no direct computation:

* **Ghost In / Ghost Out**: Act as entry and exit barriers for the density phase.
* **Stars Prep Ghosts (1-4)**: Ensure that all neighbour interactions in a given Prep loop are finished globally before the next loop begins.
* **Hydro Prep Ghost 1**: Synchronizes gas particle states if they were modified during a ``prep`` loop.


How to activate the neighbour loops
-----------------------------------

TODO (define the ``EXTRA_STARS_LOOP`` preprocessor directives at compile time)
