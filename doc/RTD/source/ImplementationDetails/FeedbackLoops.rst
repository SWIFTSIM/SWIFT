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

Extra Star Loops Activation
---------------------------

To support complex sub-grid models, the stellar feedback module supports up to four distinct preparation interaction loops. These loops enable multi-stage physics where properties are calculated and communicated across MPI ranks before the final feedback is applied to the gas.

Activating the Loops
++++++++++++++++++++

The number of interaction loops is controlled by pre-processor macros defined in the ``src/feedback.h`` header for each specific scheme. These loops must be activated to match the specific requirements of the sub-grid physics.

To activate the loops, define the following macros:

.. code-block:: c

   #define EXTRA_STAR_LOOPS_1
   #define EXTRA_STAR_LOOPS_2
   #define EXTRA_STAR_LOOPS_3
   #define EXTRA_STAR_LOOPS_4

Implementation Constraints
++++++++++++++++++++++++++

When configuring a scheme, the following rules must be observed:

* **Entry Point**: The first extra loop in the sequence must always be either ``EXTRA_STAR_LOOPS_1`` or ``EXTRA_STAR_LOOPS_2``.
* **Sequential Activation**: To maintain a coherent dependency chain in the task scheduler and MPI communication, the loops are defined in sequential order: ``prep1``, ``prep2``, ``prep3`` and ``prep4``. In the current implementation, this order cannot be changed.

Examples of Use Cases
+++++++++++++++++++++

Different models utilise these loops to handle their specific feedback mechanisms:

* **EAGLE Kinetic**: This scheme defines the first two loops (``EXTRA_STAR_LOOPS_1`` and ``EXTRA_STAR_LOOPS_2``).
* **GEAR Mechanical Feedback**:
    * **Mode 1**: Activates ``EXTRA_STAR_LOOPS_2`` and ``EXTRA_STAR_LOOPS_3``.
    * **Mode 2**: Activates ``EXTRA_STAR_LOOPS_2``, ``EXTRA_STAR_LOOPS_3``, and ``EXTRA_STAR_LOOPS_4``.
