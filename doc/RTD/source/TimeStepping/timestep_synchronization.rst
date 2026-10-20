.. Time-step synchronization
   Matthieu Schaller 9th November 2019

.. _time_step_sync:
   
Time-step synchronization
=========================

Enabled with command-line option ``--sync``.

We also implement a synchronization step to change the time-step of particles
that have been directly affected by external source terms, typically feedback
events. Durier & Dalla Vecchia (2012) showed that the Saitoh & Makino (2009)
mechanism was not sufficient in scenarios where particles receive energy in the
middle of their regular time-step. When particles are affected by feedback (see
Sections 8.1, 8.2, and 8.3 of `Schaller et al. (2024)
<https://ui.adsabs.harvard.edu/abs/2024MNRAS.530.2378S/abstract>`_), we flag
them for synchronization.

A final pass over the particles, implemented as a task acting on any cell
drifted to the current time, takes these flagged particles, interrupts their
current step to terminate it at the present time, and forces them back onto the
timeline (Section 2.4, ibid). They immediately recompute their time-step, are
assigned to an active bin for the upcoming step, and get integrated forward as
if they were on a short time-step all along. This guarantees a correct
propagation of energy and hence an efficient implementation of feedback. The
use of this mechanism is strongly recommended in simulations with external
source terms.

How it is implemented
---------------------

**Flagging.** A particle is flagged by ``timestep_sync_part()``, which sets
``limiter_data.to_be_synchronized``. It is called from the interaction
functions of the models that inject energy into a particle in the middle of its
step, for example the feedback and the black hole loops. The cell of a flagged
particle is marked (``cell_activate_sync_part()``) and the ``timestep_sync``
task of its super cell runs. The task runs after the ``timestep`` task and the
time-step limiter of the cell, and before ``kick1`` (see the chain of tasks in
:ref:`kick_drift_kick`).

**Processing.** For every flagged particle of the cell, ``runner_do_sync()``
does the following. A particle that is active at the current time has nothing
to synchronise and is only unflagged. For the others:

1. ``timestep_process_sync_part()`` cuts the current step at the current time.
   It undoes the first half-kick of the old step, then applies a kick over
   :math:`[t_{\rm beg},\, t_{\rm now}]`. The particle has now been integrated
   up to the current time, as if its step had ended here.
2. The function computes the new time-step of the particle, with the same rules
   as for any particle (see :ref:`time_step_criteria`). A pending request of
   the limiter is applied. The new bin is finally limited to
   ``e->max_active_bin``, so that the particle is guaranteed to be active in
   this step.
3. The bin and the times of the cell are updated. The ``kick1`` task then
   applies the first half-kick of the new step, as it does for any particle
   that starts a step.

The bin of a particle during the synchronisation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At the start of the process, ``timestep_process_sync_part()`` stores
``- e->min_active_bin`` in ``p->time_bin``. This value is negative and does not
represent a valid time-bin; it merely marks the particle as temporarily ready
to compute a new step. A few statements later, ``runner_do_sync()`` overwrites
it with the new bin.  Until then, any reader inspecting ``time_bin`` will see
this negative value.

This has consequences for code that reads time-bins concurrently:

* **Release builds:** A negative bin successfully passes ``part_is_active()``
  because the test is simply ``time_bin <= max_active_bin``, which correctly
  identifies a particle undergoing synchronisation as active.
* **Debug builds:** ``part_is_active()`` also performs timeline boundary
  checks. The function ``get_integer_time_end()`` returns ``0`` for a
  non-positive bin, causing debug checks to abort the run with the error:
  `"particle in an impossible time-zone! p->ti_end=0"`.

Any code that inspects time-bins simultaneously with the ``timestep_sync``
task—or reads snapshots of bins taken at that moment—must explicitly compare
``time_bin`` against ``max_active_bin`` directly (after filtering out
``time_bin_inhibited``) rather than invoking ``part_is_active()``.

Over MPI, the cells that a rank keeps for neighbouring particles are filled
with the time-bins sent by neighbours for the time-step limiter (see
:ref:`time_step_limiter`). This exchange copies every particle's bin when the
send is ready (``runner_do_pack_limiter()``) without locking the cell, waiting
only for the sender's ``timestep`` task. Because it is uncoordinated with the
``timestep_sync`` task, a foreign cell may occasionally capture the temporary
negative value. For details on how foreign cells manage this data, refer to the
developer documentation on MPI communications.
