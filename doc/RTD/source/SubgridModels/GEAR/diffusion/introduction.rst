.. GEAR FVPM metal diffusion
   Documentation added 7 August 2026

.. _gear_fvpm_diffusion_introduction:

Introduction
------------

As explained on the :ref:`gear_chemistry` page, SPH particles do not
exchange mass with one another. Feedback events (e.g. supernovae) deposit
metals onto nearby gas particles, but once deposited, those metals are
locked into the particle that received them: there is no physical mixing
process moving metal mass from one particle to its neighbours. The
:ref:`gear_smoothed_metallicity` scheme works around this by *smoothing*
the metallicity field with the SPH kernel for the purpose of, e.g.,
feedback yield tables and cooling. Crucially, it does **not** move any
metal mass between particles.

FVPM (Finite Volume Particle Method) metal diffusion is a different:
it treats each tracked chemical element as a
conserved scalar (a metal mass) that is exchanged between
neighbouring gas particles through fluxes computed on FVPM faces, and solved
with an approximate Riemann solver. The method is meshless: no explicit
mesh or tessellation is constructed, and each particle's effective volume
comes from the SPH kernel, exactly as in Gizmo's MFV/MFM hydro solvers,
which this mirrors. Here the same flux-exchange
machinery is applied only to the metal-mass field, so it works on top of
the ordinary SPH hydro (whose particles keep fixed masses) as well as on
top of Gizmo MFM hydro. The result is an explicit, conservative
diffusion/mixing model for metals, with its own diffusion coefficient,
rather than a kernel-smoothing convenience field.

Two variants of the diffusion equation are implemented: ``parabolic`` and
``hyperbolic`` equations.

Parabolic diffusion
~~~~~~~~~~~~~~~~~~~

The parabolic variant solves the standard (Fickian) diffusion equation
for the metal-mass density :math:`\rho Z`:

.. math::

   \frac{\partial(\rho Z)}{\partial t} = -\nabla\cdot\mathbf{F},
   \qquad \mathbf{F} = -K \nabla q ,

where :math:`\mathbf{F}` is the diffusive metal flux, and the
diffusivity tensor :math:`K` and the driver quantity :math:`q` are both
set by the diffusion mode in use (see `Diffusion modes`_ below). Like
all parabolic PDEs, this equation has an infinite propagation speed: an
infinitesimal perturbation at one point instantaneously affects the
solution everywhere else, which is not physical but is usually an
acceptable approximation when the diffusion coefficient and time steps
are modest compared to the resolution.

Hyperbolic diffusion
~~~~~~~~~~~~~~~~~~~~~

The hyperbolic variant instead promotes the flux :math:`\mathbf{F}` to
its own dynamical variable, coupling the same mass conservation law to a
Cattaneo (telegraph) relaxation equation:

.. math::

   \frac{\partial(\rho Z)}{\partial t} + \nabla\cdot\mathbf{F} = 0,
   \qquad
   \frac{\partial\mathbf{F}}{\partial t} = -\frac{\mathbf{F}}{\tau} - \frac{K}{\tau}\nabla q,

which introduces a finite relaxation time :math:`\tau` between the flux
and the gradient that sources it. This turns the diffusion equation into
a damped wave equation with a *finite* propagation speed
:math:`c_{\rm hyp} = \sqrt{\kappa/\tau}`, avoiding the infinite-speed
unphysicality of Fickian diffusion. As :math:`\tau \to 0`, the relaxation
equation forces :math:`\mathbf{F}\to-K\nabla q` and the system above
collapses onto the parabolic equation. The relaxation time :math:`\tau`
is itself set by a choice of mode (constant, or scaling with the local
turbulence), described together with the other runtime parameters on the
:ref:`gear_fvpm_diffusion_parameters` page.

.. note::

   Three separate mechanisms limit the diffusive flux :math:`\mathbf{F}`
   above, each guarding against a different failure mode. A **slope
   limiter** keeps the gradient reconstruction free of new local
   extrema; the code ships several classic choices (minmod, MC,
   superbee, van Leer, Koren) alongside the GIZMO-style limiter used by
   default, but switching between them is a source-code edit
   (comment/uncomment the relevant call in
   ``chemistry_slope_limiters_face.h``), not a runtime parameter. A
   **flux limiter** then acts per interaction, damping noise and
   protecting sink-like or pristine particles; its ``flux_limiter_*``
   parameters are on the :ref:`gear_fvpm_diffusion_parameters` page.
   Finally, a global **flux-corrected-transport (FCT) positivity
   limiter** rescales each particle's outgoing flux so its metal mass
   stays non-negative to within an empirical round-off margin, not an
   exact guarantee; see :ref:`gear_fvpm_diffusion_usage`.

Diffusion modes
~~~~~~~~~~~~~~~~

Both variants derive their diffusivity from one of three diffusion
modes, selected with ``diffusion_mode`` (see
:ref:`gear_fvpm_diffusion_parameters`); the mode also determines which
quantity is actually diffused, the driver :math:`q` used above:

* ``Isotropic_constant``: a constant, isotropic diffusivity
  :math:`K = \kappa\,\mathbb{1}`, with :math:`\kappa = C` a free
  parameter set directly by the user (an ordinary diffusivity, in
  length\ :sup:`2`/time). The driver is the metal *mass density*,
  :math:`q = \rho Z`.
* ``Smagorinsky``: an isotropic, turbulence-driven diffusivity
  :math:`K = \kappa\,\mathbb{1}`, with :math:`\kappa = C\,\gamma^2 h^2
  \bar\rho\,\lVert S \rVert`, with :math:`h` the particle's smoothing
  length, :math:`\gamma` the kernel's support-radius factor, so
  :math:`\gamma h` is proportional to the kernel's compact-support
  radius, and :math:`\bar\rho` the kernel-filtered density (see below).
  :math:`\kappa` grows with the local resolution, density, and the norm
  of the local shear tensor :math:`S`, i.e. stronger local
  shear/turbulence diffuses metals faster. The driver is the metal mass
  *fraction*, :math:`q = Z`.
* ``Gradient``: an anisotropic diffusion model (`Rennehan et al. 2021
  <https://doi.org/10.1093/mnras/stab1813>`_) using the
  full shear tensor rather than only its norm, :math:`K = C\,\gamma^2 h^2
  \bar\rho\,S_+`, where :math:`S_+` is :math:`S` regularised to keep only its
  positive eigenvalues (`Balarac et al. 2013
  <https://doi.org/10.1063/1.4813812>`_): a diffusivity tensor built
  from a shear tensor with negative eigenvalues would locally *sharpen*
  gradients instead of smoothing them, which this regularisation
  prevents. The driver is also the metal mass fraction, :math:`q = Z`.

In all three cases, :math:`S` is the symmetric, trace-free part of the
particle's kernel-filtered velocity gradient, with the Hubble flow added
back to its diagonal for cosmological runs before the trace is removed.
So only genuine shear/turbulence drives the diffusivity, not uniform
expansion or compression. :math:`\bar\rho` is likewise the
kernel-filtered density, not the particle's own local density.

A number of the constants entering this construction (e.g.
``GEAR_FVPM_DIFFUSION_FILTERING_SMOOTHING_FACTOR``, the weight applied
while accumulating the velocity field feeding into :math:`S`) are fixed
at compile time rather than exposed as runtime parameters; see
:ref:`gear_fvpm_diffusion_usage` for the full list.

Timestep criterion
~~~~~~~~~~~~~~~~~~~

Both variants share a chemistry CFL-like timestep parameter,
``C_CFL_chemistry`` (see :ref:`gear_fvpm_diffusion_parameters`), but the
propagation-speed difference between the two equations means they obey
genuinely different stability bounds:

* the parabolic (infinite-speed) variant is an explicit update of a
  diffusion equation, stable only for ``C_CFL_chemistry`` :math:`< 0.5`;
  its timestep also scales as :math:`(\Delta x)^2` with the local
  particle size :math:`\Delta x`, so it shrinks rapidly as the
  resolution increases or the diffusion coefficient grows;
* the hyperbolic (finite-speed) variant instead obeys an ordinary
  wave-speed CFL condition, stable for ``C_CFL_chemistry`` :math:`< 1`,
  with a timestep that scales only as :math:`\Delta x` (linear, not
  quadratic, in the resolution).

Here :math:`\Delta x` is the particle's own effective size (the FVPM/Gizmo
length scale), not the SPH smoothing length
:math:`h` used in the diffusion-mode formulas above: the two are related
resolution length scales, but not the same quantity.

.. note::

   For the full mathematical treatment (governing equations, Riemann solver, slope
   limiters, positivity limiter, timestep criteria) see the theory notes,
   which are the authoritative reference for this scheme:

   * ``theory/ParabolicDiffusion/GEAR_FVPM_Diffusion_Theory.tex``, covering both the
     parabolic and hyperbolic formulations.

Compiling the model
--------------------

The scheme is chosen at **configure time**, like every other physics
module (see the top-level project notes on the plug-in architecture), by
replacing the smoothed-metallicity scheme's ``--with-chemistry=GEAR_N``
with one of:

.. code:: bash

   --with-chemistry=GEAR-FVPM-DIFFUSION_N             # parabolic variant
   --with-chemistry=GEAR-FVPM-HYPERBOLIC-DIFFUSION_N  # hyperbolic variant

where ``N`` is the number of chemical elements to follow, same convention
as ``GEAR_N``. Both variants also work coupled to the Gizmo MFM hydro
solver instead of the default SPH hydro, by adding:

.. code:: bash

   --with-hydro=gizmo-mfm --with-riemann-solver=hllc

to the configure line. As with any other configure-time physics choice,
switching between the smoothed-metallicity scheme, the parabolic variant
and the hyperbolic variant requires reconfiguring and rebuilding; there is
no runtime switch.

The scheme is enabled automatically once the code is compiled with one of
the flags above: no extra runtime flag is required beyond what GEAR
chemistry already needs (e.g. ``--feedback`` if running with GEAR
feedback). All model parameters live under the ``GEARChemistry`` section
of the parameter file, in exactly the same section used by the smoothed
metallicity scheme; see :ref:`gear_fvpm_diffusion_parameters` for the full
list.
