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

MPI support
~~~~~~~~~~~~

.. warning::
   Running with more than one MPI rank is **not currently supported** for
   either the parabolic or the hyperbolic FVPM diffusion scheme, and this
   is enforced by an explicit, unconditional error rather than being
   merely untested.

   ``engine_maketasks()`` (:file:`src/engine_maketasks.c`) contains, for
   both ``CHEMISTRY_GEAR_FVPM_DIFFUSION`` and
   ``CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION`` builds:

   .. code:: c

      if (e->nr_nodes > 1)
        error(
            "The GEAR FVPM chemistry FCT positivity limiter does not support MPI "
            "runs yet.");

   with the code comment explaining the reason: the flux-corrected-transport
   (FCT) positivity limiter's donor-side outflow sums and per-particle
   :math:`\theta` factors are not exchanged over MPI, so a foreign
   (proxy) copy of a particle on another rank would apply a different,
   inconsistent limiter than its owning rank. The check runs
   unconditionally at every task-graph (re)build, so any run launched
   with more than one rank aborts immediately with this error; it is not
   a matter of the run silently producing wrong results with many ranks.

   Separately, the pairwise flux-exchange kernel
   (``runner_iact_chemistry_flux_exchange_decide`` in
   :file:`src/chemistry/GEAR_FVPM_DIFFUSION/chemistry_iact.h`) does contain
   a dedicated cross-rank code path that deterministically decides how to
   apply a flux when the neighbour is a foreign/proxy particle. This
   indicates the low-level interaction kernel has been written with
   multi-rank correctness in mind, but it is unreachable in practice
   because the ``engine_maketasks()`` guard above prevents any run with
   ``nr_nodes > 1`` from getting that far.

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
