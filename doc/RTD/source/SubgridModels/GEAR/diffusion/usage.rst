.. GEAR FVPM metal diffusion
   Documentation added 7 August 2026

.. _gear_fvpm_diffusion_usage:

Usage notes
-----------

Debug build: ``SWIFT_CHEMISTRY_DEBUG_CHECKS``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Configuring with

.. code:: bash

   CFLAGS="-DSWIFT_CHEMISTRY_DEBUG_CHECKS" ./configure --with-chemistry=GEAR-FVPM-DIFFUSION_N [other options]

enables extra bookkeeping and sanity checks that are too expensive to run
by default:

* Each gas particle accumulates its cumulative diffused metal mass (the
  running sum of every flux applied to it) and its cumulative
  feedback-received metal mass, per tracked element, over the whole run.
  Summing these across all particles lets you check metal-mass
  conservation directly from the snapshots (keep in mind that star
  formation *removes* metal mass from the gas budget when a gas particle
  is turned into a star).
* Two extra snapshot fields are written for gas particles,
  ``DiffusedMetalMasses`` and ``FeedbackMetalMasses`` (and, for the
  hyperbolic variant only, a third field, ``NormDiffusionFluxes``); see
  the "Gas particles" table on the :doc:`../output` page for their exact
  units and layout.
* Persistent unphysical states are reported to the run's log: if a
  particle's metal mass stays negative for many consecutive steps, or the
  total metal mass of a particle drifts to an unphysical value, a warning
  is printed (rate-limited so it does not flood the log).

These checks add overhead and are intended for validating a new setup or
debugging a suspected mass-conservation problem, not for production runs.

Fixed constants
~~~~~~~~~~~~~~~~

A handful of behaviours are set by compile-time constants in
``chemistry_properties.h``/``chemistry.h`` rather than by a runtime
parameter; changing them means editing the source and recompiling. The
ones that affect the physics or numerics a user might notice:

* ``GEAR_FVPM_DIFFUSION_FILTERING_SMOOTHING_FACTOR`` (:math:`0.8`): applied
  to the neighbour-summed part of the kernel-filtered velocity field
  feeding into the shear tensor :math:`S` (see
  :ref:`gear_fvpm_diffusion_introduction`), before the particle's own
  contribution is added, to avoid double-counting; the filtered density
  :math:`\bar\rho` used in the same formulas is not affected by this
  factor;
* ``GIZMO_SLOPE_LIMITER_BETA_MIN``/``_MAX`` (:math:`1.0`/:math:`2.0`):
  clamp the active (GIZMO-style) cell slope limiter's :math:`\beta`
  factor; the lower bound keeps the slope from being entirely zeroed out
  in smooth regions, the upper bound caps how far the reconstruction may
  extrapolate;
* ``GEAR_FVPM_DIFF_CELL_LIMITER_SHOOT_TOLERANGE`` (:math:`0.2`): the
  fractional overshoot the cell-wide positivity limiter allows, loose
  enough that it does not stall the solution;
* ``GEAR_FVPM_DIFF_WAVESPEED_ESTIMATE_DIFFERENCE_TOLERANCE``
  (:math:`10^{-8}`): floors the Riemann solver's wavespeed-difference
  denominator purely to avoid division by zero, with no physical meaning
  of its own;
* ``GEAR_FVPM_DIFF_NEGATIVITY_COUNTER_PRINT_LIMIT`` (:math:`200`): how
  many consecutive steps a particle's metal mass may stay negative
  before a console warning is triggered (production-scale logging
  cadence, not a physics threshold);
* ``GEAR_FVPM_DIFF_FLUX_LIMITER_VERBOSITY`` (off, i.e. :math:`0`, by
  default): can be turned on for extra per-interaction diagnostic
  printouts.

Two further, commented-out ``#define`` toggles in
``chemistry_properties.h`` exist purely for development use:
``GEAR_FVPM_DIFF_DEBUG_FORCE_LOOP_ONESIDED_UPDATE`` forces one-sided
pair updates for serial/thread symmetry testing (not MPI-compatible),
and ``GEAR_FVPM_DIFF_DEBUG_PAIR_VISIT_COUNT`` counts how often a
mixed-band interaction's tie-break logic takes each branch, printed at
exit.

MPI support
~~~~~~~~~~~~

.. warning::
   Running with more than one MPI rank is **not currently supported** for
   either the parabolic or the hyperbolic FVPM diffusion scheme, and this
   is enforced by an explicit, unconditional error.

   In short: **configure and run single-node/single-rank** (e.g. with
   ``--disable-mpi``, or by launching ``swift`` rather than ``swift_mpi``,
   or ``swift_mpi`` on a single rank) when using either FVPM diffusion
   scheme.

Worked examples
~~~~~~~~~~~~~~~~

A set of standalone examples exercising this scheme live under
``examples/ChemistryTests/``. Each has its own ``README`` and ``run.sh``;
consult those for the exact configure line and any run-specific caveats.

* ``MetalDiffusionOnePeak``: a non-cosmological, hydro-only homogeneous
  box in which a single gas particle at the box centre is seeded with
  non-zero metallicity while every other particle starts metal-free.
  Exercises constant isotropic diffusion of a point-like source. Its
  ``README`` documents important caveats for the hyperbolic variant at
  large :math:`\tau`: in that wave-dominated regime the exact solution is
  a thin expanding shell rather than a smooth blob, so visible front
  smearing at finite resolution is an expected outcome of the test, not a
  defect.
* ``MetalDiffusionTwoPeaks``: the same homogeneous box, but with two
  metal-seeded gas particles placed symmetrically at :math:`x=L/4` and
  :math:`x=3L/4`. Exercises diffusion from two simultaneous point
  sources.
* ``MetalDiffusionDiscontinuity``: a box with a sharp step in the initial
  metallicity (one half of the box at one metal mass fraction, the other
  half at a different one), rather than a point source. Runnable in 1D,
  2D or 3D (set via ``--with-hydro-dimension`` and the ``dimension``
  environment variable to ``run.sh``); a good test of how the scheme
  resolves a sharp gradient rather than a near-delta seed.
* ``MetalDiffusionWaveInterference``: reuses ``MetalDiffusionTwoPeaks``'
  initial conditions, but because :math:`L/4` and :math:`3L/4` are exactly
  antipodal on the periodic box, the two seeds' finite-speed fronts meet
  at two simultaneous, symmetric crossing points under the **hyperbolic**
  scheme. This is a direct demonstration that hyperbolic diffusion behaves
  as a genuine finite-speed wave equation; it is not meaningful for the
  parabolic scheme, which has no front to cross.
* ``MetalAdvectionTestGEAR``: tests that metals are correctly advected
  with the mass fluxes exchanged by mass-flux hydro schemes (e.g. Gizmo
  MFV), for the GEAR chemistry element set. Configured with
  ``--with-hydro=gizmo-mfv --with-riemann-solver=hllc``; to see a visible
  advection effect, the hydro solver needs to be run in Eulerian mode
  (``GIZMO_FIX_PARTICLES`` in ``src/const.h``), since in Lagrangian mode
  little mass is exchanged in the first place.
* ``MetalAdvectionTestEAGLE``: the same advection test, for the EAGLE
  chemistry scheme instead of GEAR.
