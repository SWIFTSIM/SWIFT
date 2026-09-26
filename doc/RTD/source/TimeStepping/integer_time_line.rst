.. Integer time-line
   Matthieu Schaller 9th November 2019

.. _integer_time_line:

Integer time-line & time-bins
=============================

SWIFT does not keep track of time with floating-point numbers. Every
quantity that decides *when* something happens (the start and the end of a
time-step, whether a particle is active, when to write a snapshot) is an
integer on a fixed grid, the *integer time-line*. This page describes the
grid, the *time-bins* that give a time-step to every particle and the
conversion to physical time. The definitions are in ``src/timeline.h``.

The time-line
-------------

The integer time is stored in a 64-bit integer (``integertime_t``). The
whole simulation, from the start to the end time, is divided into
:math:`2^{57}` equal intervals. Each one is called a *tick*.

.. note::

   **A tick is the unit of the integer time.** It is the smallest step that
   time can make in the code: nothing happens between two ticks. The integer
   time is a number of ticks counted from the start of the run (it is 0 at the
   start and :math:`2^{57}` at the end). A time-step is a whole number of
   ticks. The physical length of a tick is called ``time_base`` (its value is
   given below). For example, a run that covers 10 Gyr has ticks of about 2
   seconds.

The total number of ticks is the constant ``max_nr_timesteps``, which is
``1 << (num_time_bins + 1)`` with ``num_time_bins = 56``.

For a non-cosmological run, the length of one tick is

.. math::

   \Delta t_{\rm tick} = \frac{t_{\rm end} - t_{\rm begin}}{2^{57}}

and this value is stored in ``time_base`` (``engine_init()``). The physical
time of the integer time :math:`t_{\rm int}` is
:math:`t_{\rm begin} + t_{\rm int}\,\Delta t_{\rm tick}` (see
``engine_step()``).

For a cosmological run, the variable that is divided into ticks is the
logarithm of the scale-factor. The value of ``time_base`` is
:math:`(\log a_{\rm end} - \log a_{\rm begin}) / 2^{57}` (``cosmology_init()``
in ``cosmology.c``). The physical time and the scale-factor of every step are
computed from the integer time by ``cosmology_update()``. As a consequence,
``dt_min`` and ``dt_max`` are then intervals of :math:`\log a` (see
:ref:`Parameters_time_integration`). The time-step criteria of the particles
return physical time-steps. They are multiplied by the Hubble rate
:math:`H(a)`, which is stored in ``cosmology->time_step_factor``, before they
are converted to ticks.

The simulation is finished when the integer time reaches :math:`2^{57}`
(``engine_is_done()``).

Time-bins
---------

Each particle carries a *time-bin*, the field ``time_bin``. A particle in the
time-bin :math:`b \geq 1` uses time-steps that are

.. math::

   \Delta t_{\rm int}(b) = 2^{b+1}\ {\rm ticks}
   = \left(t_{\rm end} - t_{\rm begin}\right) \times 2^{\,b-56}

long. This is the function ``get_integer_timestep()``. It returns 0 for a
bin :math:`b \leq 0`, which is not a valid step length. The shortest possible
time-step is the one of the bin 1 and has 4 ticks. The longest one is the one
of the bin 56 and is the whole time-line. The function ``get_time_bin()`` goes
the other way. It returns :math:`\lfloor\log_2 \Delta t_{\rm int}\rfloor - 1`,
so a time-step that is not a power of two is *rounded down* to the bin below.
The ``+ 1`` in the exponent of the length exists to keep the *half* of the
shortest step an integer. The kick-drift-kick integrator (see
:ref:`kick_drift_kick`) needs the middle of every step.

.. figure:: figures/time_bins.png
   :width: 100%
   :alt: Length of the time-step of every time-bin

   The length of the time-step of the 56 time-bins, in ticks (left) and for
   the Sedov blast example (right). In this example the run is :math:`0.05`
   long, ``dt_min`` is :math:`10^{-9}` and ``dt_max`` is :math:`10^{-2}`. Every
   bin doubles the length of the step. SWIFT prints the two circled values at
   the start of a run, as "Minimal timestep size (on time-line)" and "Maximal
   timestep size (on time-line)". They are the lengths
   :math:`(t_{\rm end}-t_{\rm begin})/2^k` that are just below ``dt_min`` and
   ``dt_max``.

A few other values of the time-bin mark special particles. They are not valid
step lengths.

* ``time_bin_inhibited`` (``num_time_bins + 2``): the particle has been
  removed from the simulation (``part_is_inhibited()``).
* ``time_bin_not_created`` (``num_time_bins + 3``): the slot is a spare one,
  kept for a particle that will be created during the run.
* ``time_bin_not_awake`` (``-num_time_bins``): this is not a time-bin. It is
  the value of ``limiter_data.wakeup`` when no neighbour asked for the
  particle to be woken up (see :ref:`time_step_limiter`).
* A negative ``time_bin`` is also written for a short time while a particle
  is synchronised on the time-line (see :ref:`time_step_sync`).

Which bins are active
---------------------

The steps of the bin :math:`b` start and end at the integer times that are
multiples of :math:`\Delta t_{\rm int}(b)`. Every length is a power of two.
Hence a multiple of a long step is also a multiple of all the shorter ones. At
the end of a step of the bin :math:`b`, all the particles in the bins
:math:`1,\ldots,b` are also ending a step. This is why SWIFT never needs a
list of the active particles. It only needs the highest active bin at the
current time. A particle is active when

.. code-block:: c

   p->time_bin <= e->max_active_bin

which is the whole content of ``part_is_active()``. The highest active bin is
computed from the integer time by ``get_max_active_bin()``. It is the position
of the lowest set bit of the integer time, minus one. The lowest active bin is
the bin whose step is as long as the step that has just been taken
(``get_min_active_bin()``).

.. figure:: figures/active_bins.png
   :width: 100%
   :alt: Active time-bins as a function of the integer time

   The bins that end a step at each integer time (blue), for the eight
   shortest bins, and the resulting highest active bin (bottom). The bin 1 is
   active at every point. The bin :math:`b` is active once every
   :math:`2^{b-1}` of the points. The highest active bin follows the sequence
   1, 2, 1, 3, 1, 2, 1, 4, and so on.

The time-line is shared by all the particles. Two particles are on the same
step only when the integer time is a multiple of both of their step lengths.
The choice of the new bin of a particle is therefore not free. A particle can
only take a longer step if its current time is a multiple of that longer step
(see :ref:`time_step_criteria`).
