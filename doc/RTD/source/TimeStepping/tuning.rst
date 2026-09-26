.. Reading the table of the steps and tuning the time-steps

.. _time_step_tuning:

Reading the table of the steps
==============================

SWIFT prints one line per step, on the standard output and in the file
``timesteps.txt``. This page explains the columns, shows what a real table
looks like, and uses it to explain what the time-step parameters cost. The
example is the ``SedovBlast_3D`` example
(``examples/HydroTests/SedovBlast_3D``), with its default parameter file, the
:math:`64^3` glass and four threads. The table that is used for the figures is
kept in ``TimeStepping/figures/sedov_3d_timesteps.txt`` together with the
scripts that draw them.

The columns
-----------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Column
     - Meaning
   * - ``Step``
     - The number of the step. Step 0 is the initialisation, at the start time,
       where all the time-bins are active.
   * - ``Time``
     - The time at the end of the step, in internal units. In a cosmological run
       it is the cosmic time.
   * - ``Scale-factor``, ``Redshift``
     - The values at the end of the step (1 and 0 without cosmology).
   * - ``Time-step``
     - The time between the previous step and this one (``e->time_step``). It is
       the length of the step of the lowest active bin. It is *not* the step of
       every particle.
   * - ``Time-bins``
     - Two numbers: the lowest and the highest active time-bin of the step
       (see :ref:`integer_time_line`).
   * - ``Updates``, ``g-Updates``, ``s-Updates``, ``Sink-Updates``, ``b-Updates``
     - The number of active gas, gravity, star, sink and black hole particles in
       the step, summed over all the ranks. These are the particles that took
       a kick2 and a new time-step.
   * - ``Wall-clock time``
     - The duration of the step, in milliseconds.
   * - ``Props``
     - The properties of the step, as the sum of the flags that apply:
       1 tree rebuild, 2 redistribution of the particles, 4 repartitioning of
       the domain, 8 statistics, 16 snapshot, 32 restart file, 64 structure
       finding, 128 friends-of-friends, 256 mesh gravity, 512 power spectra. For
       instance 25 is a rebuild, statistics and a snapshot.
   * - ``Dead time``
     - The average, over the threads (and ranks), of the time that a thread
       spent without a task to run during the step, in milliseconds.

How to read it
--------------

.. figure:: figures/sedov_steps.png
   :width: 100%
   :alt: Length of the step, fraction of updated particles and active bins
         of the Sedov blast run

   The steps of the Sedov blast run. **Top:** the time between steps. It grows
   from :math:`2.4\times 10^{-5}` to :math:`3.9\times 10^{-4}` as the blast
   spreads and the particles slow down. **Middle:** the fraction of the
   particles that are updated. Most steps update a small fraction. The curves are
   the steps that end the same bin: a longer bin activates more
   particles. The diamonds are the *synchronised* steps, where all the
   particles are active. **Bottom:** the lowest and highest active bins. The
   highest active bin sets how many particles are updated: the peaks are the
   steps where a long bin ends.

The run has 386 steps after the initialisation and 8 synchronised steps, and
uses the bins 45 to 56. Some numbers that the table gives:

* **The gain of individual time-steps.** The shortest step is
  :math:`2.44\times 10^{-5}` and the run is :math:`0.05` long. A code with a
  single global step of that length would need 2048 steps and 2048 updates of
  every particle. In the run, the sum of the updates is 19.2 times the number
  of particles, that is about 107 times fewer.
* **A step costs something even with few active particles.** 238 of the 386
  steps update less than 1% of the particles, and together they take 11% of the
  wall-clock time. The 146 steps that update fewer than 100 particles, when they
  have no rebuild and no output, take 3.2 ms each (median).
* **The dead time is large in these steps.** The median ratio of dead time to
  wall-clock time is 33% for the steps that update fewer than 100 particles,
  and 2% for the steps that update more than :math:`10^4`. A step with few
  active particles cannot keep all the threads busy.
* **Rebuilds and outputs are expensive.** The 33 steps with a rebuild (``Props``
  contains 1) are 18% of the wall-clock time. A step with a rebuild and fewer
  than 1000 updates costs about 60 ms (median). The 4 steps that write a snapshot
  cost about 0.9 s.
* **The full steps.** The 8 synchronised steps update all the particles and take
  24% of the time.

.. figure:: figures/sedov_step_cost.png
   :width: 75%
   :alt: Wall-clock time of a step against the number of updated particles

   The cost of the steps of the same run. Above about :math:`10^3` updated
   particles, the cost grows roughly in proportion to the number of updates
   (about 6 ms per :math:`10^3` particles for the regular steps of more than
   :math:`10^4` updates). Below that, the cost is a floor of a few
   milliseconds. Steps with a
   tree rebuild (squares) and steps that write a snapshot (diamonds) are far
   above the trend.

What the parameters do to the cost
----------------------------------

``dt_max``
   It caps the longest time-step. If it is much shorter than what the particles
   ask for, every particle is forced to a short step, the bins collapse
   towards the same few values, and the run takes many more steps with many
   active particles. SWIFT rounds it down to a length of the form
   :math:`(t_{\rm end} - t_{\rm begin})/2^k` and prints the result as "Maximal
   timestep size (on time-line)".

``dt_min``
   It does not make steps shorter or longer. It is a safeguard: a particle that
   asks for less stops the run (see below). Set it low enough not to stop a
   healthy run, but not so low that a runaway particle can drive the run
   through millions of tiny steps before somebody notices.

The spread of the bins
   The cost of a run is the number of steps times the fixed cost of a step,
   plus the cost of the updates. The shortest bin that is populated sets the
   number of steps. A few particles in a very short bin (a blast, a very dense
   region, a sink) make the run take steps in which only they are active. Those
   steps cost the fixed part and use few threads (see the dead time above). The
   time-step limiter and the synchronisation add to it, since they shorten the
   steps of the neighbours of such particles (see :ref:`time_step_limiter`).

The number of outputs and rebuilds
   Every output step drifts all the cells and writes the data (see
   :ref:`kick_drift_kick`). Every rebuild costs more than a regular step. For
   runs with self-gravity, the parameter ``Gravity:rebuild_active_fraction``
   (see :ref:`Parameters_gravity`) can trigger a rebuild when the fraction
   of active gravity particles is large.

The particle time-step criteria (``SPH:CFL_condition``, ``Gravity:eta``,
``SPH:max_volume_change`` and the parameters of the subgrid models) decide the
step that every particle asks for (see :ref:`time_step_criteria`). They are the
main way to change the accuracy and the cost together.

Messages and errors
-------------------

At the start of a run, SWIFT prints the length of the tick and the limits that
it found on the time-line:

.. code-block:: none

   engine_config: Absolute minimal timestep size: 3.469447e-19
   engine_config: Minimal timestep size (on time-line): 7.450581e-10
   engine_config: Maximal timestep size (on time-line): 6.250000e-03

The first value is ``time_base``, the length of a tick (see
:ref:`integer_time_line`). The two other values are explained in the caption of
the figure in that page. SWIFT stops at the start with an error when:

* ``dt_min`` is larger than ``dt_max``: "Minimal time-step size (...) must be
  smaller than maximal time-step size (...)";
* ``dt_min`` is smaller than the length of a tick: "Minimal time-step size
  smaller than the absolute possible minimum dt=...";
* in a non-cosmological run, ``dt_max`` is larger than the length of the run:
  "Maximal time-step size larger than the simulation run time t=...".

During a run:

* "part (id=...) wants a time-step (...) below dt_min (...)" (and the same for
  gpart, spart, bpart and sink) means that the criteria of one particle ask for
  a step below ``dt_min``. This is the safeguard doing its job. Look at the
  state of the particle with the identifier that is printed (its velocity,
  sound speed, acceleration, smoothing length and softening). The cause is
  often a physical problem rather than a wrong parameter.
* "Not limiting particle with id ... because it needs to be synced" is a warning
  of the limiter (see :ref:`time_step_limiter`). It only means that a particle
  was to be woken up by a neighbour and to be synchronised in the same step, and
  that the synchronisation takes care of it.
