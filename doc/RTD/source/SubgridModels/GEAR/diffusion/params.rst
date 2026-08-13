.. GEAR FVPM metal diffusion
   Documentation added 7 August 2026

.. _gear_fvpm_diffusion_parameters:

Model parameters
-----------------

All parameters of the FVPM diffusion scheme live in the same
``GEARChemistry`` section of the parameter file as the smoothed
metallicity scheme. Most parameters are **shared** between the parabolic
and hyperbolic variants; a small set applies to the **hyperbolic variant
only** (marked explicitly below), because it does not exist as a physical
quantity in the parabolic equation.

Initial metallicity
~~~~~~~~~~~~~~~~~~~~

The first two parameters are identical to the smoothed-metallicity
scheme and are not repeated in full here; see
:ref:`gear_smoothed_metallicity`:

* ``initial_metallicity``: the initial (non-smoothed) mass fraction of
  each tracked element for all particles. If negative, the initial
  metallicities are instead read from the initial conditions. This
  parameter is mandatory.
* ``scale_initial_metallicity`` (optional, default ``0``): whether to
  scale the initial metallicity of each element using the solar
  composition from the feedback yield tables.

Diffusion coefficient and mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``diffusion_coefficient`` (optional, default :math:`1.0`): the
  normalisation constant :math:`C` that sets the diffusion coefficient
  :math:`\kappa`, i.e. :math:`\kappa \propto C`. Its physical meaning
  depends on ``diffusion_mode`` below: for ``Isotropic_constant`` it *is*
  :math:`\kappa` directly (an ordinary diffusivity, length\ :sup:`2`/time);
  for ``Smagorinsky`` and ``Gradient`` it is a dimensionless prefactor
  multiplying a local turbulence/shear estimate, giving :math:`\kappa` a
  mass-diffusivity form instead.
* ``diffusion_mode`` (mandatory string): selects how :math:`\kappa` (or
  the diffusion matrix :math:`K`) is computed. Must be one of:

  * ``Isotropic_constant``: constant, isotropic diffusion with coefficient
    :math:`\kappa = C`.
  * ``Smagorinsky``: isotropic, Smagorinsky-type turbulent diffusion,
    proportional to the local shear-tensor norm :math:`|S|`.
  * ``Gradient``: anisotropic diffusion following the gradient model of
    `Rennehan et al. (2021) <https://doi.org/10.1093/mnras/stab1813>`_,
    proportional to the shear tensor :math:`S` itself rather than only
    its norm.

  Any other value raises a configuration error at start-up.

.. warning::

   :math:`\kappa`'s units depend on ``diffusion_mode`` (see above), but
   the ``DiffusionMatrices`` snapshot output (see
   :ref:`gear_fvpm_diffusion_usage`) always reports the *converted*
   ordinary diffusivity :math:`D`, not the raw internal matrix, so its
   units are the same regardless of which mode is in use.

HLL Riemann solver parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These control the HLL-type approximate Riemann solver used to compute the
metal-mass flux across each FVPM face (`Hopkins 2017
<https://arxiv.org/abs/1602.07703>`_ conventions):

* ``hll_riemann_solver_psi`` (optional, default :math:`0.1`): the
  :math:`\psi` parameter controlling the solver's wavespeed estimate. Must
  be non-negative; a negative value raises a configuration error.
* ``hll_riemann_solver_epsilon`` (optional, default :math:`0.0`): a
  tolerance parameter for the solver. Must satisfy
  :math:`0 \leq \varepsilon \leq 1`; values outside this range raise a
  configuration error.

Flux limiter parameters
~~~~~~~~~~~~~~~~~~~~~~~~

These control the per-interaction flux limiter applied on top of the raw
Riemann-solver flux, before the global FCT positivity limiter:

* ``flux_limiter_noise_gate`` (optional, default :math:`10^{-15}`):
  fluxes below this fraction of the source particle's metal mass are
  treated as numerical noise and zeroed, which prevents a slow
  accumulation ("ratchet effect") of spurious tiny fluxes.
* ``flux_limiter_safety`` (optional, default :math:`0.5`): a global safety
  factor applied by the rational flux limiter.
* ``flux_limiter_sink_stability`` (optional, default :math:`0.25`): caps
  the fractional change a single neighbour interaction may cause to a
  particle's metal mass in one step, protecting sink-like particles from
  large single-interaction swings.
* ``flux_limiter_startup`` (optional, default :math:`0.1`): the maximum
  fraction of the source particle's mass that may flow into a pristine
  (:math:`Z=0`) neighbour in one step, which seeds diffusion into
  metal-free regions without an unbounded initial jump.

Timestep parameter
~~~~~~~~~~~~~~~~~~~

* ``C_CFL_chemistry`` (mandatory): the Courant-Friedrich-Levy-like
  constant used to compute the chemistry diffusion timestep constraint.
  This parameter is required in both the parabolic and hyperbolic
  variants, but its stable range depends on which one is compiled: the
  **parabolic** variant requires ``C_CFL_chemistry`` :math:`< 0.5`,
  while the **hyperbolic** variant is stable up to
  ``C_CFL_chemistry`` :math:`< 1`. Values at or above these bounds make
  the run unstable; see :ref:`gear_fvpm_diffusion_introduction` for why
  the two bounds differ.

Hyperbolic-only parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following parameters exist **only** when compiled with
``--with-chemistry=GEAR-FVPM-HYPERBOLIC-DIFFUSION_N``; they configure the
relaxation time :math:`\tau` that turns the parabolic equation into the
Cattaneo/telegraph equation:

* ``tau`` (mandatory): the physical relaxation time. For
  ``relaxation_time_mode: Constant`` this *is* :math:`\tau`; for
  ``relaxation_time_mode: Turbulent`` it is instead a multiplicative
  factor applied to a turbulence-based estimate of :math:`\tau`. **It must
  not** be set to :math:`0.0`: a zero relaxation time defeats the purpose
  of hyperbolic diffusion (it would reduce to the parabolic limit while
  still paying the extra hyperbolic bookkeeping) and raises a
  configuration error.
* ``relaxation_time_mode`` (mandatory string): selects how :math:`\tau` is
  computed. Must be one of:

  * ``Constant``: :math:`\tau` is the fixed value given by ``tau``.
  * ``Turbulent``: :math:`\tau` scales with the inverse of the local
    turbulent shear, :math:`\tau \propto 1/\|S\|`, with ``tau`` acting as
    the proportionality factor.

  Any other value raises a configuration error.
* ``hyperbolic_limiter_scope`` (optional string, default
  ``AllComponents``): which components of the HLL solver's dissipative
  flux the causal bound rescales. Must be one of:

  * ``Density``: rescale only the mass-density (monopole) dissipation
    term.
  * ``AllComponents``: also rescale the dissipation of the flux-vector
    components, not only the density term.

  Any other value raises a configuration error.

Example parameter sections
~~~~~~~~~~~~~~~~~~~~~~~~~~~

A minimal ``GEARChemistry`` section for the **parabolic** variant:

.. code:: YAML

   GEARChemistry:
     initial_metallicity: 0.0           # Initial metal mass fraction (or < 0 to read from the ICs)
     scale_initial_metallicity: 0       # (Optional) Scale initial metallicity by the solar composition. (Default: 0)
     diffusion_mode: Isotropic_constant # Isotropic_constant, Smagorinsky or Gradient
     diffusion_coefficient: 1.0         # (Optional) Normalisation constant for kappa. (Default: 1.0)
     hll_riemann_solver_psi: 0.1        # (Optional) HLL solver wavespeed parameter, >= 0. (Default: 0.1)
     hll_riemann_solver_epsilon: 0.0    # (Optional) HLL solver tolerance, in [0, 1]. (Default: 0.0)
     flux_limiter_noise_gate: 1e-15     # (Optional) Noise-gate threshold. (Default: 1e-15)
     flux_limiter_safety: 0.5           # (Optional) Global flux-limiter safety factor. (Default: 0.5)
     flux_limiter_sink_stability: 0.25  # (Optional) Sink-particle stability cap. (Default: 0.25)
     flux_limiter_startup: 0.1          # (Optional) Startup fraction into pristine neighbours. (Default: 0.1)
     C_CFL_chemistry: 0.4               # CFL-like constant for the chemistry diffusion timestep.
                                        # Must be < 0.5 (parabolic) or < 1 (hyperbolic), see above.

The **hyperbolic** variant additionally requires ``relaxation_time_mode``
and ``tau``, and accepts the optional ``hyperbolic_limiter_scope``:

.. code:: YAML

   GEARChemistry:
     # ... all the parabolic parameters above, plus:
     relaxation_time_mode: Constant          # Constant or Turbulent
     tau: 0.1                                # Relaxation time (Constant mode) or its prefactor (Turbulent mode). Must not be 0.
     hyperbolic_limiter_scope: AllComponents # (Optional) Density or AllComponents. (Default: AllComponents)
