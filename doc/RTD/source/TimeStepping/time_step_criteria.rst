.. How the time-step of a particle is chosen
   Darwin Roduit, 2026

.. _time_step_criteria:

Choosing the time-step of a particle
====================================

When a particle ends a step, the ``timestep`` task gives it a new time-bin
(see :ref:`kick_drift_kick`). This page explains how the new time-step is
computed. There are two stages. First, every physics module that acts on the
particle proposes a *physical* time-step, and the smallest proposal is kept.
Second, this time-step is converted to the integer time-line and adjusted, so
that it fits the rules of the time-line (see :ref:`integer_time_line`).

Stage 1: the proposals
----------------------

The proposals depend on the type of the particle. The functions that combine
them are in ``src/timestep.h``.

Gas particles (``get_part_timestep()``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The new time-step is the minimum of:

* The hydrodynamical time-step, ``hydro_compute_timestep()``. It is set by the
  hydro scheme (in ``src/hydro/<scheme>/hydro.h``). It is typically a CFL
  condition. For the SPHENIX scheme it is
  :math:`\Delta t = 2\,\gamma_{\rm k}\,{\rm CFL}\; a\, h / (a_{\rm cs}\, v_{\rm sig})`,
  where :math:`\gamma_{\rm k}` is the ratio of the kernel support to :math:`h`,
  :math:`v_{\rm sig}` is the signal velocity of the particle, :math:`a` and
  :math:`a_{\rm cs}` are cosmological factors (equal to one without
  cosmology) and the CFL factor is the parameter ``SPH:CFL_condition``.
* The magneto-hydrodynamical time-step, ``mhd_compute_timestep()``.
* The cooling time-step, ``cooling_timestep()``, if cooling is switched on.
* The gravity time-steps of the particle, if it has a gravity counterpart:
  ``gravity_compute_timestep_self()`` and ``external_gravity_timestep()``. The
  first one is the criterion of Gadget-2 (type 0),
  :math:`\Delta t = \sqrt{2\,\eta\, a\,\epsilon/|{\bf g}|}`, where
  :math:`{\bf g}` is the physical acceleration, :math:`\epsilon` the softening
  length, converted to its Plummer-equivalent value, and :math:`\eta` is the
  parameter ``Gravity:eta``.
* The time-step of the forcing terms, ``forcing_terms_timestep()``.
* The time-step of the chemistry model, ``chemistry_timestep()`` (for instance
  metal diffusion).

Two further limits are applied to that minimum.

* The change of the smoothing length: the time-step cannot let :math:`h` change
  by more than a factor set by the parameter ``SPH:max_volume_change`` (the
  default is :math:`1.4`). The limit is :math:`|\log(f_h)\,h/\dot h|`, where
  :math:`f_h` is that maximal factor for :math:`h` (``log_max_h_change``).
* In cosmological runs only, the constraint on the displacement of the
  particles, ``dt_max_RMS_displacement``. It is computed by
  ``engine_recompute_displacement_constraint()``, and is set by the parameters
  ``TimeIntegration:max_dt_RMS_factor`` and ``TimeIntegration:dt_RMS_use_gas_only``.

In a cosmological run, the result is then multiplied by the Hubble rate
(``cosmology->time_step_factor``), so that it is an interval of :math:`\log(a)`.
Finally it is limited by ``dt_max``. If the result is below ``dt_min``, SWIFT
stops with the error

.. code::

   "part (id=...) wants a time-step (...) below dt_min (...)".

This protects against a run that could never finish.

Other particles
~~~~~~~~~~~~~~~

Gravity-only particles (``get_gpart_timestep()``) use only the gravity
time-steps. The constraint on the displacement and the conversion for the
cosmology are the same as for gas. Star, black hole and sink particles add
the time-step of their own model to the gravity ones:

* Stars (``get_spart_timestep()``): ``stars_compute_timestep()``, plus the
  radiative transfer time-step if it is used. The parameters that limit the
  time-step of young and old stars are described in :ref:`Parameters_Stars`.
* Black holes (``get_bpart_timestep()``): ``black_holes_compute_timestep()``.
* Sinks (``get_sink_timestep()``): ``sink_compute_timestep()``.

Each of them stops with an error if the result is below ``dt_min``.

A gas particle that has a gravity counterpart gives its time-bin to it, so
that both always use the same steps.

Stage 2: the rules of the time-line
-----------------------------------

The function ``make_integer_timestep()`` turns the physical time-step into a
valid number of ticks. It applies these rules in this order.

1. **Conversion and rounding.** The time-step is converted to ticks and rounded down to a time-bin (``get_time_bin()``).

2. **Neighbours.** For gas particles, the bin cannot be more than ``time_bin_neighbour_max_delta_bin`` (that is, 2 by default) bins above the smallest time-bin of the neighbours of the particle. In other words, a particle cannot have a step more than four times longer than any of its neighbours. The smallest bin of the neighbours is collected in the density loop (``runner_iact_timebin()``, field ``limiter_data.min_ngb_time_bin``). Other particle types have no such limit.

3. **Growth.** If the old bin is positive, the new step cannot be longer than twice the old step. A step can thus grow by one bin at a time, but it can shrink by any amount.

4. **Alignment.** A step of length :math:`T` can only be *longer* than the previous one if the current time :math:`t` is an absolute multiple of :math:`T` from the start of the simulation (:math:`t=0`). If it is not, the particle keeps its previous step for one more step. (A shorter step always fits, because every multiple of a long step is a multiple of the short ones.)

When the radiative transfer with sub-cycling is used, the step of a gas
particle is finally limited to ``max_nr_rt_subcycles`` times its radiative
transfer step (see the pages on radiative transfer).

Worked Example
~~~~~~~~~~~~~~

To see how these rules interact sequentially, consider a particle currently on **bin 3 (16 ticks)** whose physical criteria request **100 ticks**, with the smallest neighbour bin set to **2**. 

The table below traces how this request is progressively constrained depending on whether the step happens to end at :math:`t = 48` or :math:`t = 64`.

.. list-table::
   :header-rows: 1
   :widths: 35 32 33

   * - Constraint Rule
     - Evaluation at :math:`t = 48`
     - Evaluation at :math:`t = 64`
   * - **1. Request** (100 ticks)
     - **Bin 5** (64 ticks)
     - **Bin 5** (64 ticks)
   * - **2. Neighbours** (:math:`\le \text{min\_ngb} + 2`)
     - Capped at :math:`2 + 2` :math:`\rightarrow` **Bin 4** (32 ticks)
     - Capped at :math:`2 + 2` :math:`\rightarrow` **Bin 4** (32 ticks)
   * - **3. Growth** (:math:`\le 2 \times \text{old bin}`)
     - Allowed (32 :math:`\le 2 \times 16`) :math:`\rightarrow` **Bin 4**
     - Allowed (32 :math:`\le 2 \times 16`) :math:`\rightarrow` **Bin 4**
   * - **4. Alignment** (:math:`T` must divide :math:`t`)
     - 32 does not divide 48 :math:`\rightarrow` **Reverts to Bin 3** (16 ticks)
     - 32 divides 64 cleanly :math:`\rightarrow` **Accepts Bin 4** (32 ticks)

As this example demonstrates, particles do not instantly jump to their requested criteria. They are strictly bounded by their neighbours, growth limits, and global timeline alignment, ensuring system-wide synchronization.

The time-step limiter and the synchronisation (see :ref:`time_step_limiter`
and :ref:`time_step_sync`) can shorten the step of a particle in the middle of
its step, whatever the rules above have decided.
