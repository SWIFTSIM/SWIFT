.. Time-step synchronization
   Matthieu Schaller 9th November 2019

.. _time_step_sync:
   
Time-step synchronization
=========================

Enabled with command-line option ``--sync``.

We also implement a synchronization step to change the time-step 
of particles that have been directly affected by external source 
terms, typically feedback events. Durier & Dalla Vecchia (2012) 
showed that the Saitoh & Makino (2009) mechanism was not sufficient 
in scenarios where particles receive energy in the middle of their 
regular time-step. When particles are affected by feedback (see 
Sections 8.1, 8.2, and 8.3 of Schaller et al. (2024), MNRAS 530:2), 
we flag them for synchronization. A final pass over the particles, 
implemented as a task acting on any cell which was drifted to the 
current time, takes these flagged particles, interrupts their current 
step to terminate it at the current time and forces them back onto 
the timeline (Section 2.4, ibid) at the current step. They then recompute 
their time-step and get integrated forward in time as if they were 
on a short time-step all along. This guarantees a correct propagation 
of energy and hence an efficient implementation of feedback. The use 
of this mechanism is always recommended in simulations with external source terms.

How it is implemented
---------------------

**Flagging.** A particle is flagged by ``timestep_sync_part()``, which sets
``limiter_data.to_be_synchronized``. It is called from the interaction
functions of the models that inject energy in a particle in the middle of its
step, for example the feedback and the black hole loops. The cell of a flagged
particle is marked (``cell_activate_sync_part()``) and the ``timestep_sync``
task of its super cell runs. The task runs after the ``timestep`` task and the
time-step limiter of the cell, and before ``kick1`` (see the chain of tasks
in :ref:`kick_drift_kick`).

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
   ``e->max_active_bin``, so that the particle is active in this step.
3. The bin and the times of the cell are updated. The ``kick1`` task then
   applies the first half-kick of the new step, as it does for any particle
   that starts a step.

The bin of a particle during the synchronisation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At the end of the first step, ``timestep_process_sync_part()`` stores
``-e->min_active_bin`` in ``p->time_bin``. The value is negative and is not
the bin of a valid step. It only marks the particle as ready to compute a new
step. A few statements later, ``runner_do_sync()`` replaces it by the new bin.
Until then, any reader of ``time_bin`` sees the negative value.

This has a consequence for the code that reads time-bins. In a release build,
a negative bin passes ``part_is_active()``, since ``part_is_active()`` is the
test ``time_bin <= max_active_bin``, and that is the desired answer for a
particle that is being synchronised. In a debug build, ``part_is_active()``
also checks the end of the step of the particle. The function
``get_integer_time_end()`` returns 0 for a bin that is not positive, and the
check stops the run with the error "particle in an impossible time-zone!
p->ti_end=0". Any code that can read a bin at the same time as the
``timestep_sync`` task of that cell, or that reads a copy of the bins taken at
such a time, must therefore compare ``time_bin`` with ``max_active_bin``
directly (after excluding the ``time_bin_inhibited`` value) instead of calling
``part_is_active()``.

Over MPI, the cells that a rank keeps for the particles of its neighbours are
filled with the time-bins that the neighbours send, for the time-step limiter
(see :ref:`time_step_limiter`). The exchange copies the bin of every particle of
the cell when the send is ready (``runner_do_pack_limiter()``). The copy is
done without a lock on the cell, and the only task of the sender that it
waits for is the ``timestep`` task. It is not ordered with respect to the
``timestep_sync`` task. A foreign cell can therefore hold the negative value
above for some of its particles. The page on the MPI communication in the
developer documentation describes the foreign cells and the data that they
hold.
