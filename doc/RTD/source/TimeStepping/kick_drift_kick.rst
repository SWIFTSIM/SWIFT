.. The time integration scheme and the time-line of a step

.. _kick_drift_kick:

The integration scheme and the steps of the code
================================================

SWIFT integrates the particles with a kick-drift-kick (leapfrog) scheme in
which every particle has its own time-step (see :ref:`integer_time_line` and
:ref:`time_step_criteria`). This page describes what the scheme does to a
particle, which tasks do it, and how the code decides when the next step
happens.

The scheme for one particle
---------------------------

Consider a particle that starts a step at the integer time
:math:`t_{\rm beg}` and takes a time-step :math:`\Delta t` (a whole number of
ticks, given by its time-bin). It ends the step at
:math:`t_{\rm end} = t_{\rm beg} + \Delta t`. The scheme has three operations.

* **Kick** (``kick_part()`` and its counterparts for the other particle
  types). It changes the velocity, and other quantities, using the current
  acceleration. A step has two kicks, each of half the step. The first kick
  (*kick1*) covers :math:`[t_{\rm beg},\, t_{\rm beg} + \Delta t/2]` and the
  second (*kick2*) covers :math:`[t_{\rm beg} + \Delta t/2,\, t_{\rm end}]`.
  The kicked velocity is stored in ``xpart->v_full``.
* **Drift** (``drift_part()``). It moves the particle to a given time, by
  :math:`x \leftarrow x + v_{\rm full}\,\Delta t_{\rm drift}`. It also predicts
  the velocity ``v`` and the other fields that the neighbour loops use at that
  time. The drift is not tied to the steps of the particle. It goes from the
  last drift of the cell to the current time (see below).
* **Force**. The neighbour loops compute the accelerations at the end of the
  step.

The order of the operations for a particle is: *kick1* at the start of the
step, the drifts to the times at which the particle is needed, the force at
the end of the step, then *kick2*. The acceleration computed at
:math:`t_{\rm end}` is used twice: for the *kick2* of the step that ends and for
the *kick1* of the next step. This is why, at the end of a step, the
particle first takes its new time-step, and only then its *kick1* (see the
tasks below).

Without cosmology, the interval of a kick or of a drift is the difference of the
two integer times multiplied by ``time_base`` (``kick_get_grav_kick_dt()`` and
the other ``kick_get_*_kick_dt()`` functions). In a cosmological run these
functions return instead the cosmological factors between the two integer times
(``cosmology_get_grav_kick_factor()``, ``cosmology_get_drift_factor()`` and the
equivalent functions for the hydro and thermal kicks), so that the scheme is
the same in comoving coordinates.

.. figure:: figures/kdk_timeline.png
   :width: 100%
   :alt: Kick-drift-kick timeline of a particle over two steps

   A particle takes a step of the bin 3 (16 ticks) and then a step of the bin 4
   (32 ticks), to scale. The growth to the bin 4 is possible because the end of
   the first step, :math:`t=32`, is a multiple of 32 (see
   :ref:`time_step_criteria`). The kicks split every step in two halves. The
   drift uses the velocity that the *kick1* has produced, which is valid at
   the middle of the step (filled circle). The *kick2* brings the velocity to
   the end of the step (open circle). At the tick where a step ends, the tasks
   run in the order: force, *kick2* of the step that ends, ``timestep`` (which
   chooses the new bin), *kick1* of the step that starts.

"Active" and "starting"
-----------------------

A particle is *active* at the integer time :math:`t` when it ends a step at
:math:`t`. It is *starting* when it begins a step at :math:`t`. At the moment
that it takes its new step, it is both. The two tests are
``part_is_active()`` and ``part_is_starting()``. They give the same answer
(``time_bin <= max_active_bin``), because the particles that end a step are
the ones that start a new one. They are used at different moments of the step
and check different things in the debug builds (the end of the step is not in
the past for "active", the start is not in the future for "starting").

* The loops over the particles that run before the new time-bin is chosen use
  "active": the neighbour loops, *kick2* and the ``timestep`` task.
* *kick1* runs after the ``timestep`` task and uses "starting" with the
  *new* time-bin.

For a cell, ``cell_is_active_hydro()`` is true when the earliest end of a step
in the cell is the current time (``hydro.ti_end_min == ti_current``), and
``cell_is_starting_hydro()`` when the latest start of a step in the cell is
the current time (``hydro.ti_beg_max == ti_current``). The other particle
types have the same functions. A particle that was removed from the
simulation (``time_bin_inhibited``) is neither.

The tasks of a step
-------------------

For each cell that holds the tasks of the particles (the *super cell*), the
following chain runs every step. Only the parts of the chain that are needed
by an active cell are activated.

.. code-block:: none

   drift -> density, gradient and force loops -> ghosts -> end_force
         -> kick2 -> timestep -> limiter loops -> timestep_limiter
         -> timestep_sync -> kick1

* ``drift`` (``runner_do_drift_part()``) is only active if a task needs the
  particles of the cell at the current time.
* ``end_force`` finishes the acceleration of the active particles.
* ``kick2`` (``runner_do_kick2()``) applies the second half-kick to the active
  particles.
* ``timestep`` (``runner_do_timestep()``) calls ``get_part_timestep()`` and the
  equivalent functions, gives every active particle its new time-bin, and
  records the times of the cell (below). Subgrid tasks such as star formation
  run between *kick2* and ``timestep``.
* The limiter loops, ``timestep_limiter`` and ``timestep_sync`` exist only with
  the options ``--limiter`` and ``--sync``. They can change the bin of
  particles (see :ref:`time_step_limiter` and :ref:`time_step_sync`).
* ``kick1`` (``runner_do_kick1()``) applies the first half-kick of the new
  step to the starting particles.
* ``timestep_collect`` runs at the level of the top cells, after the
  ``timestep`` tasks of all their super cells and, when they exist, the limiter
  and synchronisation tasks. It collects the times of the cells.

Note that the first half-kick of a step is done at the *end of the previous
step*, in the same call of ``engine_step()``. Only the *drift* and the force
belong to the following call.

Drifts are lazy
~~~~~~~~~~~~~~~

Inactive particles are not moved at every step. Each cell remembers the time
of its last drift (``hydro.ti_old_part``). When a task needs the particles of
a cell at the current time, the unskip phase flags the cell and the ``drift``
task moves all of them from ``ti_old_part`` to the current time
(``cell_drift_part()``). The cells that no task needs stay behind and are moved
later, over a longer interval. All the cells are drifted at once
(``engine_drift_all()``) when the tree is rebuilt and to write an output.

The times of the cells and the next step
----------------------------------------

Each cell stores two times per particle type: ``ti_end_min`` (the earliest
end of a step of any of its particles) and ``ti_beg_max`` (the latest start
of a step). The ``timestep`` task computes them for the super cell. For an
active particle they come from the new step. For an inactive one they come
from its current step (``get_integer_time_end()`` and
``get_integer_time_begin()``). The ``timestep_collect`` task and then
``engine_collect_end_of_step()`` combine them up to the top cells and, over
MPI, across the ranks. The result is

.. code-block:: c

   e->ti_end_min = min(ti_hydro_end_min, ti_gravity_end_min, ti_sinks_end_min,
                       ti_stars_end_min, ti_black_holes_end_min);

the time of the next step: the earliest moment at which any particle ends a
step.

One call of ``engine_step()`` takes this sequence.

1. The time moves forward: ``ti_old = ti_current``,
   ``ti_current = ti_end_min``. The highest and the lowest active bins are set
   from these two times (``get_max_active_bin()`` and
   ``get_min_active_bin()``), and the physical time and the step size
   ``e->time_step`` are updated.
2. ``engine_prepare()`` activates the tasks that are needed in this step
   (``engine_unskip()``). If the tree must be rebuilt or the domain
   repartitioned, everything is first drifted to the current time, the
   rebuild is done, and the tasks are activated after that.
3. The tasks run (``engine_launch()``). Over MPI, the times of the top-level
   cells are then exchanged (``engine_synchronize_times()``).
4. ``engine_collect_end_of_step()`` collects the counts of updated particles
   and the times of the cells, and computes the next ``ti_end_min``.
5. The outputs are written if one is due (``engine_io()``, see below).

The run stops when ``ti_current`` reaches the end of the time-line.

The first step (number 0 in the table of the steps) is special. It runs the
initialisation (``engine_init_particles()``) at :math:`t=0`, where all the bins
are active, and gives every particle its first time-step.

Outputs do not change the steps
-------------------------------

Snapshots, statistics and the other outputs happen at times that are on the
integer time-line (``ti_next_snapshot``, ``ti_next_stats`` and so on). They do
**not** shorten the time-steps of the particles. At the end of a step,
``engine_io()`` checks whether the next end of step (``ti_end_min``) is later
than the next output. If it is, the engine temporarily sets the current time to
the time of the output, marks all the bins as inactive and drifts every
particle to that time, and writes the output. Then it puts the time back. The
particles are not kicked, so their steps and bins are unchanged. The price is a
drift of all the cells and the writing itself. In the table of the steps, the
column ``Props`` shows which steps did I/O (see :ref:`time_step_tuning`).

Synchronised steps
~~~~~~~~~~~~~~~~~~

A step at which the highest active bin is at least as high as the highest bin
that contains particles is *synchronised*: all the particles are active. It
happens when the integer time is a multiple of the longest step in use. There
are few of them, at most one per longest step. Some quantities can only be
computed on such steps. For instance, the exact gravity forces that are used
to check the accuracy (see the parameter ``only_when_all_active`` in the
section on the gravity force checks of the parameter file documentation) need
all the particles to be at the same time.
