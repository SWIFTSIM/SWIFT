.. Time-step limiter
   Matthieu Schaller 9th November 2019

.. _time_step_limiter:

Time-step limiter
=================

SWIFT uses a two-part time-step limitation strategy to ensure energy
conservation (section 3.4, `Schaller et al. (2024)
<https://ui.adsabs.harvard.edu/abs/2024MNRAS.530.2378S/abstract>`_).

1. **Upcoming step restriction (Always active):** When a particle computes its
   next time-step, it restricts its size relative to its neighbors. This
   baseline mechanism requires no extra tasks and is always enabled.
2. **Active wake-up and interruption (Enabled with** ``--limiter`` **):** Based on
   Saitoh & Makino (2009), this optional extension wakes up *inactive*
   particles whose active neighbors are running on much smaller time-steps,
   interrupting their ongoing steps if necessary.

What the limiter does
---------------------

The first limit we impose is to limit the time-step of active particles.
When a particle computes the size of its next time-step, typically using the CFL condition,
it also additionally considers the time-step size of all the particles
it interacted with during the loop computing accelerations. We then demand
that the particle of interest’s time-step size is not larger than a
factor :math:`\Delta` of the minimum of all the neighbours’ values. We typically
use :math:`\Delta = 4` which fits naturally within the binary structure of the
time-steps in the code. This first mechanism is always activated in
Swift and does not require any additional loops or tasks; it is, however,
not sufficient to ensure energy conservation in all cases.

The time-step limiter proposed by Saitoh & Makino (2009) is
also implemented in SWIFT and is a recommended option for all
simulations not using a fixed time-step size for all particles (enabled via ``--limiter``). This
extends the simple mechanism described above,
by also considering inactive particles and waking them up
if one of their active neighbours uses a much smaller time-step size.

As shown by Saitoh & Makino (2009), this is necessary to conserve energy and hence
yield the correct solution even in purely hydrodynamics problems
such as a Sedov–Taylor blast wave. The additional loop over the
neighbours is implemented by duplicating the already existing tasks
and changing the content of the particle interactions to activate the
requested neighbours.

How it is implemented
---------------------

The two mechanisms use two fields of ``struct timestep_limiter_data``, which
every gas particle carries (``src/timestep_limiter_struct.h``). The factor
:math:`\Delta = 4` is two time-bins:
``time_bin_neighbour_max_delta_bin`` in ``src/timeline.h`` is 2 and a step of
the bin :math:`b + 2` is four times longer than one of the bin :math:`b`.

* ``min_ngb_time_bin`` is the smallest time-bin among the neighbours. It is
  collected in the density loop (``runner_iact_timebin()``), and it implements
  the first mechanism. When a particle chooses its new step, the step is
  limited to a bin at most two above ``min_ngb_time_bin`` (see
  :ref:`time_step_criteria`).
* ``wakeup`` is the request to wake the particle up. It is ``time_bin_not_awake``
  when nobody asked.

The wake-up requests are made by the limiter loop. This loop runs for the
active particles, after their new time-bin has been chosen (see the chain of
tasks in :ref:`kick_drift_kick`). For an active particle :math:`i` and each
neighbour :math:`j`, ``runner_iact_nonsym_limiter()`` does

.. code-block:: c

   if (pj->time_bin > pi->time_bin + time_bin_neighbour_max_delta_bin)
      accumulate_max_c(&pj->limiter_data.wakeup, -pi->time_bin);

Because the largest of the negative values wins, ``wakeup`` holds minus the
smallest bin among the active neighbours that asked. The new bin of the woken
particle is ``-wakeup + 2``: two bins above that neighbour, so that the
condition is fulfilled.

The ``timestep_limiter`` task (``runner_do_limiter()``) then applies the
requests to the particles of a cell (``timestep_limit_part()``). A particle
that was flagged for a synchronisation (see :ref:`time_step_sync`) and that is
not active is skipped, with the warning "Not limiting particle with id ...
because it needs to be synced". There are two cases.

* The particle is active. It is ending a step now, and its new bin, computed
  by the ``timestep`` task, is too long. The bin is replaced by ``-wakeup + 2``
  and the new step starts at the current time.
* The particle is inactive. It is in the middle of a longer step, which must be
  *interrupted*. This is the interesting case, and it is described below.

Interrupting a step
~~~~~~~~~~~~~~~~~~~

Let :math:`t` be the current time. The step that is interrupted goes from
:math:`t_{\rm beg}` to :math:`t_{\rm end}`. Its first half-kick has already
been applied, up to the middle of the step :math:`t_{\rm mid}`. The particle is
given the shorter step :math:`\Delta t_{\rm new}` of the new bin. The new step
starts at :math:`t_{\rm beg,new}`, which is the latest point of the form
:math:`t_{\rm beg} + k\,\Delta t_{\rm new}` that is not after :math:`t`.

To see how this works in practice, consider the concrete timeline illustrated
in the figure below: a particle on bin 5 (:math:`\Delta t = 64` ticks) is
interrupted at :math:`t = 20` by an active neighbour on bin 1, forcing a switch
to bin 3 (:math:`\Delta t_{\rm new} = 16` ticks, starting at :math:`t_{\rm
beg,new} = 16`).

The function executes three precise adjustments:

1. **Undo the old half-kick:** It reverses the first half-kick of the old step
   over :math:`[t_{\rm beg}, t_{\rm mid}]` using a negative interval. *(In the
   figure: undoing the kick over* :math:`[0, 32]` *)*.
2. **Shift to the new start:** It applies a kick over :math:`[t_{\rm beg},\,
   t_{\rm beg,new}]` to bring the velocity forward to the start of the new
   step. *(In the figure: applying the kick over* :math:`[0, 16]` *)*.
3. **Catch up on the new half-kick:** If the new bin is not active at the
   current time, it applies the missing first half-kick of the new step over
   :math:`[t_{\rm beg,new},\, t_{\rm beg,new} + \Delta t_{\rm new}/2]`. If the
   bin is active, the standard ``kick1`` task handles it instead. *(In the
   figure: applying the missing half-kick over* :math:`[16, 24]` *)*.

The particle now ends its new, shorter step correctly at :math:`t_{\rm
beg,new} + \Delta t_{\rm new}`.

.. figure:: figures/limiter_interrupt.png
   :width: 100%
   :alt: A particle interrupted by the time-step limiter

   A particle of the bin 5 (64 ticks) has done the first half of its step
   (blue, top). At :math:`t=20`, an active neighbour of the bin 1 wakes it up.
   The new bin is :math:`1 + 2 = 3` (16 ticks), and the latest start of a step
   of 16 ticks that is not after :math:`t` is 16. The function (1) undoes the
   kick over :math:`[0, 32]`, (2) applies the kick over :math:`[0, 16]` and (3)
   applies the missing first half-kick of the new step over :math:`[16, 24]`,
   because the bin 3 is not active at :math:`t=20`. The particle now ends its
   step at :math:`t=32`. The figure is drawn with the same integer arithmetic
   as ``timestep_limit_part()``.

Cost, MPI, and Usage
--------------------

The limiter needs the extra loop and extra tasks, which is why it is an option
(``--limiter``). Several of the run modes switch it on together with the
synchronisation (for instance ``--eagle``, ``--flamingo`` and ``--gear``). When
a particle is woken, its step becomes shorter than what its own criteria asked
for, so the number of active particles per step goes up. The bins of the
particles of a region with a strong event (a blast, a feedback event) are then
spread over a wide range (see :ref:`time_step_tuning`).

Over MPI, the limiter loop also needs the time-bins of the particles on the
neighbouring ranks. A separate exchange (``task_subtype_limiter``) sends them
once the ``timestep`` task of the sending cell is done. See
:ref:`time_step_sync` for what a rank can find in the bins it receives.
