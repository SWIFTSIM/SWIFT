.. MAGMA SPH
   Matthieu Schaller, November 2026

.. _magma_sph:

MAGMA2 SPH
==========

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Contents:

This is SWIFT's implementation of the SPH scheme of the MAGMA2 code
(`Rosswog 2020 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.498.4230R>`_,
MNRAS 498, 4230). It is a density-energy SPH scheme (the thermal variable is
the internal energy per unit mass) whose gradients are computed by matrix
inversion instead of kernel gradients, and whose artificial dissipation acts
on slope-limited reconstructions of the velocity and internal energy at the
pair mid-point rather than on the raw particle differences. The scheme
combines:

+ The standard kernel-sum density and SWIFT's usual smoothing-length
  iteration on the neighbour number (``SPH:resolution_eta``);
+ Linearly exact "gradient functions" :math:`\mathbf{G}_a` built from the
  correction matrix :math:`\mathsf{C}_a` (Rosswog 2020, eqs. 4-6),
  computed in a third (gradient) loop over the neighbours, with a
  regularised inverse for planar or filamentary neighbourhoods and a
  per-pair fallback to kernel gradients;
+ Three conservative formulations of the equations of motion, MI1 and MI2
  from the paper and MI3, a variant of MI2 introduced in SWIFT, selected at
  configure time (see below);
+ Artificial viscosity (:math:`\beta = 2`) acting on the difference of
  velocities reconstructed to the pair mid-point with a van Leer slope
  limiter and a short-distance cut-off (eqs. 14-23), with the time-dependent
  per-particle :math:`\alpha_a` (maximum 1, floor 0.5) and the shear limiter of
  `Cullen & Dehnen (2010) <https://ui.adsabs.harvard.edu/abs/2010MNRAS.408..669C>`_;
+ Artificial conduction (eq. 24) with the paper's velocity-based signal
  speed plus a term driven by the pressure jump left after reconstruction,
  and a per-particle coefficient in the style of SPHENIX that keeps the
  conduction off in shocks and around feedback events;
+ The Durier & Dalla Vecchia (2012) time-step limiter and the usual SWIFT
  entropy and pressure floors.

To configure with this scheme, use

.. code-block:: bash

   ./configure --with-hydro=magma2 --with-kernel=quintic-spline

``--with-hydro=magma2`` is the MI3 formulation. ``magma2-MI1``,
``magma2-MI2`` and ``magma2-MI3`` select a formulation explicitly. The
``Scheme`` attribute of the ``HydroScheme`` group of the snapshots records
which one was used (``MAGMA-2-MI1``, ``MAGMA-2-MI2`` or ``MAGMA-2-MI3``).

The three formulations
----------------------

With :math:`\mathbf{d}_{ab} = \mathbf{r}_b - \mathbf{r}_a`,
:math:`V_b = m_b / \rho_b` and the correction matrix

.. math::
   \mathsf{C}_a = \Big[\sum_b V_b\, \mathbf{d}_{ab} \otimes \mathbf{d}_{ab}\,
   W_{ab}(h_a)\Big]^{-1},

each particle owns a gradient function for every neighbour,
:math:`\mathbf{G}_a = \mathsf{C}_a \mathbf{d}_{ab} W_{ab}(h_a)` and
:math:`\mathbf{G}_b = \mathsf{C}_b \mathbf{d}_{ab} W_{ab}(h_b)`, which are
exact for linear fields and odd under the exchange of the two particles.
Writing :math:`Q_a`, :math:`Q_b` for the viscous pressures and
:math:`\mathbf{v}_{ab} = \mathbf{v}_a - \mathbf{v}_b`:

.. math::
   \text{MI1 (eqs. 2-3):}\quad
   \frac{d\mathbf{v}_a}{dt} = -\sum_b m_b \left[\frac{P_a + Q_a}{\rho_a^2}
   \mathbf{G}_a + \frac{P_b + Q_b}{\rho_b^2} \mathbf{G}_b\right], \qquad
   \frac{du_a}{dt} = \sum_b m_b \frac{P_a + Q_a}{\rho_a^2}\,
   \mathbf{v}_{ab} \cdot \mathbf{G}_a

.. math::
   \text{MI2 (eqs. 10-11):}\quad
   \frac{d\mathbf{v}_a}{dt} = -\sum_b m_b \frac{P_a + P_b + Q_a + Q_b}
   {\rho_a \rho_b} \mathbf{G}_{ab}, \qquad
   \frac{du_a}{dt} = \sum_b m_b \frac{P_a + Q_a}{\rho_a \rho_b}\,
   \mathbf{v}_{ab} \cdot \mathbf{G}_{ab}, \qquad
   \mathbf{G}_{ab} = \tfrac{1}{2}(\mathbf{G}_a + \mathbf{G}_b)

.. math::
   \text{MI3:}\quad
   \frac{d\mathbf{v}_a}{dt} = -\sum_b m_b \left[\frac{P_a + Q_a}{\rho_a \rho_b}
   \mathbf{G}_a + \frac{P_b + Q_b}{\rho_a \rho_b} \mathbf{G}_b\right], \qquad
   \frac{du_a}{dt} = \sum_b m_b \frac{P_a + Q_a}{\rho_a \rho_b}\,
   \mathbf{v}_{ab} \cdot \mathbf{G}_a

MI3 keeps the :math:`1/(\rho_a \rho_b)` weighting of MI2 (the
:math:`\sigma = 1` member of the paper's eq. 8, which is what removes the
surface tension of MI1 at contact discontinuities, see the square test in
the paper) but attaches each pressure to the gradient function of its own
particle, as MI1 does. All three conserve linear momentum exactly and the
energy equation is the exact partner of the momentum equation; angular
momentum is conserved only approximately in all of them because
:math:`\mathbf{G}_a` is not parallel to :math:`\mathbf{d}_{ab}`.

Which one to use:

+ **MI3 (default)**: MI2's freedom from surface tension with the full
  divergence operator of every particle. The motivation is the behaviour at
  large smoothing-length contrasts: when :math:`h_a \gg h_b` most neighbours
  :math:`b` do not see :math:`a` (:math:`W_{ab}(h_b) = 0`), so
  :math:`\mathbf{G}_b = 0` and MI2's average collapses to
  :math:`\mathbf{G}_a / 2`; the hot centre of a Sedov blast then cools too
  slowly in MI2. In MI3 only the :math:`P_b` term vanishes, as it should.
+ **MI2**: the paper's recommended formulation; use it to reproduce the
  paper's MI2 results. Expect a hotter blast centre than MI3 at high
  :math:`h` contrast.
+ **MI1**: the paper's Newtonian limit of the Rosswog (2015) formulation.
  It has the surface tension of standard SPH at contact discontinuities
  (the square becomes a diamond) and is provided for comparison; it also
  reaches the largest fraction of pair fallbacks because the angle test is
  applied to :math:`\mathbf{G}_a` and :math:`\mathbf{G}_b` separately.

Gradient functions and fallbacks
--------------------------------

The gradient loop accumulates :math:`\mathsf{C}_a^{-1}`, the matrix
gradients of the velocity, internal energy and pressure (eq. 18), the
Cullen & Dehnen shock-presence sum, the SPH Laplacian of :math:`u` and the
largest sound speed and viscosity coefficient among the neighbours. In the
extra ghost task the matrix is inverted with a regularisation: directions in
which the condition number exceeds ``SPH:gradient_max_condition_number``
(planar or filamentary neighbourhoods) are clipped while the well-sampled
directions keep their exact inverse. A particle falls back to standard SPH
(kernel gradients, no reconstruction) when the matrix cannot be inverted at
all or when :math:`h > 0.99\, h_{\max}`; the reason is written to the
``FallbackFlags`` snapshot field. In the force loop, a pair falls back to
kernel gradients when the gradient function points away from the pair
separation or makes an angle larger than 0.5 rad with it (MI1 and MI3 test
:math:`\mathbf{G}_a` and :math:`\mathbf{G}_b` separately, MI2 their
average); the number of such pairs is written to ``GradientFallbackPairs``.

Because of the :math:`h_{\max}` fallback, **set** ``SPH:h_max``
**generously**: particles whose smoothing length is capped interact with
plain SPH.

Reconstruction and slope limiter
--------------------------------

Every dissipative term acts on the jump left after each particle's gradient
has been extrapolated to the pair mid-point,
:math:`\tilde q_a = q_a + \tfrac{1}{2} \Phi_{ab} (\nabla q)_a \cdot
\mathbf{d}_{ab}`, with the limiter of eqs. 21-23: the van Leer factor
:math:`\min(1, 4A/(1+A)^2)` of the ratio :math:`A` of the two particles'
directional slopes (no reconstruction when the slopes have opposite signs),
times an exponential cut-off that switches the reconstruction off for pairs
closer than ``SPH:limiter_eta_crit`` smoothing lengths (default: the mean
inter-particle separation, :math:`1/\eta`) with a width of
:math:`0.1\,\gamma_K` smoothing lengths, :math:`\gamma_K` being the ratio of
the kernel support to :math:`h`; these are the paper's values translated to
SWIFT's :math:`h` convention (the paper's :math:`h` is half the kernel
support). The cut-off is scaled by the convergence fraction of the flow,
:math:`(\nabla \cdot \mathbf{v})^2 / ((\nabla \cdot \mathbf{v})^2 +
|\mathsf{S}|^2)` with :math:`\mathsf{S}` the traceless strain (the larger
of the two particles' values for a pair), so that close
pairs in pure shear (differential rotation) keep their reconstruction and
the viscosity does not act on the shear. The velocity, internal energy and
pressure each have their own limiter; the internal-energy reconstruction is
additionally replaced by the pair mean if it would invert the sign of the
jump.

The reconstruction is linear: it uses the first-order term of the paper's
eq. 17 only. The quadratic term (eqs. 19-20) is not implemented, by design.

Artificial viscosity
--------------------

The viscous pressure is the paper's eqs. 14-15, :math:`Q_a = \rho_a(-\alpha
c_a \mu_a + \beta \mu_a^2)` with :math:`\mu_a = \min(0, \tilde{\mathbf{v}}_{ab}
\cdot \boldsymbol{\eta}_a / (\eta_a^2 + \epsilon^2))` built from the
reconstructed velocity difference, :math:`\boldsymbol{\eta}_a =
\mathbf{r}_{ab} / h_a` and :math:`\epsilon = 0.1`. :math:`\beta` and
:math:`\epsilon` are compile-time constants at the paper's values;
:math:`\alpha` is the particle's time-dependent coefficient :math:`\alpha_a`
described below, with the paper's :math:`\alpha = 1` as its maximum.
:math:`Q` is only applied to pairs that approach (raw velocities,
including the Hubble flow) and only where its heating term
:math:`\mathbf{v}_{ab} \cdot \mathbf{G}` is positive. Note that
:math:`\mu_a` uses SWIFT's :math:`h`, see the departures below.

Each particle's :math:`Q_a` is multiplied by its Cullen & Dehnen (2010)
shear limiter :math:`\xi_a` (their eqs. 17-18): the convergence fraction
with the divergence boosted by the shock-presence factor :math:`2(1-R_a)^4`,
where :math:`R_a = \rho_a^{-1} \sum_b \mathrm{sign}(\nabla\cdot\mathbf{v}_b)
m_b W_{ab}` (the sum includes :math:`a` itself) is :math:`-1` when every
neighbour converges and :math:`+1` when every neighbour diverges. In
disordered flow with mixed signs :math:`R_a` scatters about zero by its
particle noise, the particle's own term :math:`m_a W(0)/\rho_a` plus about
:math:`1/\sqrt{N_{\rm ngb}}` (0.44 with the quintic spline at
:math:`\eta = 1.5` in 2D, 0.21 in 3D). The factor therefore reads
:math:`R_a` through a dead zone of 1.1 times that noise (capped at 0.9),
computed per particle from its mass, density and smoothing length:
:math:`|R_a|` below
it is mapped to 0 (the factor 2) and the rest rescaled so that
:math:`R_a = \mp 1` keep the factors 32 and 0. Without the dead zone the
converging half of the particles of a cold Keplerian disc get a factor of up
to 6 and the limiter lets them plunge; with it the disc holds and the
factor keeps its full strength in a front. The limiter is floored at
``SPH:viscosity_limiter_floor`` (default 0) and acts on the linear term of
:math:`Q` only for pairs with :math:`|\mu_a| > 0.5\, c_a`: subsonic shear
layers keep the full linear viscosity that maintains their particle order,
cold supersonic shear (Keplerian discs) is relieved of it.

The coefficient :math:`\alpha_a` is the time-dependent one of Cullen &
Dehnen (2010) in the form used by SPHENIX: the source :math:`A_a = \xi_a H^2
\max(0, -d(\nabla\cdot\mathbf{v})/dt)` of a converging flow sets
:math:`\alpha_{\rm loc} = \alpha_{\max} A_a/(A_a + c_s^2)`, which is of order
:math:`(kH)^2` times the amplitude in a linear wave and 1 in a shock;
:math:`\alpha_a` jumps up to :math:`\alpha_{\rm loc}` and decays back
towards it exponentially over :math:`H / (2\, \ell\, c_s)` with
:math:`\ell =` ``SPH:viscosity_length`` (default 0.05, ten sound crossings of
the support), and is floored at ``SPH:viscosity_alpha_min`` (default 0.5).
:math:`\alpha_{\max} = 1` is a compile-time constant. Particles start at
the floor. The same coefficient is the shock measure of the conduction cap
below. Note that the limiter :math:`\xi_a` acts twice, in the trigger and on
:math:`Q`.

The floor is the viscosity of smooth flow. The reconstruction leaves a few
per cent of the velocity jump of a pair to :math:`Q` (gradient noise and the
first-order mid-point reconstruction), so the scheme is not inviscid in
smooth flow and the residual is damped at the floor's strength. The floor
was set in 2026-11 on the scheme-comparison suite: against the paper's
constant :math:`\alpha = 1`, a floor of 0.3 halves the damping of a
linear sound wave (4.5 to 2 per cent per period at 32 particles per
wavelength; the damping is linear in the floor) and the floor of 0.5 cuts
it to 3 per cent, the floor cuts the numerical viscosity of the
Taylor-Green vortex (by 2.7 at 0.3, by 1.7 at 0.5), improves the Gresho
vortex and the growth rate of the Kelvin-Helmholtz instability, and doubles
the extent of the inertial range of driven subsonic turbulence at
:math:`64^3` (the :math:`k = 50` to 100 band triples and reaches 80 per cent
of Arepo's at the same particle number). Shocks, the Zeldovich pancakes,
the Keplerian rings and Galilean invariance are unchanged. The price is
particle noise at contact discontinuities, which the constant
:math:`\alpha` was damping, and the floor is the lever on it: the density
error of the square test is 0.06 at constant :math:`\alpha`, 0.09 at a
floor of 0.5 and 0.12 at 0.3 (SPHENIX 0.08), the dense mass leaking out of
the square 1, 6 and 7 per cent, and the Evrard collapse and the blob test
follow the same order. A floor of 0.1 adds visible noise in shear layers.
The default of 0.5 was chosen after the 2026-11 conduction scan, which
found the conduction constants at their optimum and the floor the only
remaining lever on the contacts: it trades the last 40 per cent of the
smooth-flow gain for half of the contact noise and a third of the Evrard
loss. ``SPH:viscosity_alpha_min: 1`` restores the paper's constant
viscosity in the force loop; the conduction cap (below) reads the
coefficient before the floor and is unchanged.

Artificial conduction
---------------------

The conduction term is the paper's eq. 24 with the signal speed

.. math::
   v^{\rm sig}_{ab} = \alpha_u\, |\tilde{\mathbf{v}}_a - \tilde{\mathbf{v}}_b|
   + \alpha_P \sqrt{\frac{2\, \max(0, |\Delta P_{\rm rec}| - f_0 \bar P)}
   {\rho_a + \rho_b}},

where the first term is the paper's gravity-safe choice (eq. 26,
:math:`\alpha_u = 0.05`) and the second (:math:`\alpha_P = 1`,
:math:`f_0 = 0.05`, all three compile-time constants) acts on the pressure
jump that the slope-limited reconstruction of the thermal pressure does
*not* remove: nothing across a resolved pressure gradient (hydrostatic
equilibrium), the full kernel-scale pressure blip across an entropy jump or
in a shock. The thermal pressure :math:`\rho u (\gamma - 1)` is used here
even when a pressure floor is active, so that the floor does not drive
conduction at its boundary. The pair term is scaled by the pressure-weighted
mean of the two particles' coefficients :math:`\alpha^D_a`
(``DiffusionParameters`` in the snapshots), which grow at ``SPH:diffusion_beta``
:math:`H \nabla^2 u / \sqrt{u}` where :math:`u` has a local minimum, decay
towards ``SPH:diffusion_alpha_min`` over a signal crossing time of the
support, and are capped by :math:`1 - (\alpha^{\rm ngb}_{\max} -
\alpha_{\min}) / (\alpha_{\max} - \alpha_{\min})`, the SPHENIX limiter, with
:math:`\alpha^{\rm ngb}_{\max}` the largest viscosity coefficient of the
particle (its current value) and its neighbours. Here :math:`\alpha_{\min}`
is the smaller of ``SPH:viscosity_alpha_min`` and its default of 0.5, and
the coefficient is read before the floor of the force loop is applied, so
that raising the floor does not weaken the cap: at a floor of 1 the force
loop has the paper's constant viscosity while the cap still switches the
conduction off in shocks and around heated particles, as at the default.
With the default floor of 1 the coefficient is this cap alone: the paper's
full conduction wherever the gas is not viscous, none in shocks and around a
particle that a feedback event has just heated (the feedback hook raises its
:math:`\alpha_a` to the maximum). Setting ``SPH:diffusion_alpha_min: 0``
gives SPHENIX's behaviour, conduction only where the Laplacian source has
switched it on. The decay's signal speed is the sum of the particle's and
its hottest neighbour's sound speeds, both without the pressure floor, as
the conduction terms themselves.

Time step
---------

Each particle takes the minimum of the acceleration criterion (eq. 35,
:math:`\sqrt{2H/|\mathbf{a}|}`) and the Courant criterion (eq. 36) in the
form

.. math::
   \Delta t = C_{\rm CFL}\, \frac{H_a}{c^{\max}_a + 0.6\,\alpha_{\max}\,(c^{\max}_a
   + 2 \tilde\mu_a)},

with :math:`H_a` the kernel support, :math:`c^{\max}_a` the largest sound
speed of the particle and its neighbours (the paper uses a global time step
set by the hottest particle), :math:`\alpha_{\max} = 1` the bound of the
per-particle viscosity coefficient (which reaches it within a step of a
shock forming, before the particle's next time-step computation) and
:math:`\tilde\mu_a` the largest softened approach velocity over the
particle's pairs. SWIFT's criteria are written
with the kernel support :math:`H`, the paper's with :math:`h = H/2`, so
**the recommended** ``SPH:CFL_condition: 0.2`` **is the paper's** :math:`C =
0.4` and SWIFT's generic 0.1 is the paper's 0.2.

The choice of the criterion and of the CFL number was tested in 2026-11 on
16 tests of the scheme-comparison suite (shocks, smooth flows, cold flows,
interfaces, self-gravity) with a scan over :math:`C_{\rm CFL} = 0.05` to
0.4, for the form above, for the same form with the particle's own sound
speed in place of :math:`c^{\max}_a`, for the pair form of the other SWIFT
schemes, :math:`\max_b (c_a + c_b - \beta \mu_{ab})`, and for SPHENIX, with
and without the time-step limiter. The outcome:

+ With the form above the errors are set by the resolution, not by the
  time step, up to :math:`C_{\rm CFL} = 0.2` on every test; the Noh implosion
  is the most demanding one and its density error is 8 per cent larger at
  0.2 than at 0.1. At 0.3 the strong blast wave, the shearing shock and Noh
  degrade by several to tens of per cent and the smooth 2D tests sit at their
  stability edge; at 0.4 particles collapse onto each other. **Use 0.2;**
  0.3 only for 3D shock problems where a 30 per cent saving in steps is
  worth a few per cent of accuracy.
+ The two cheaper forms of the signal speed take 2.5 to 3 times fewer steps
  at the same :math:`C_{\rm CFL}` but degrade on every shock from 0.2 and fail
  from 0.3. At equal accuracy they need the same number of steps as the form
  above (their 0.1 is its 0.2 on the Sod and blast tests), except next to a
  hot region, where they are cheaper and lose energy conservation: on the
  Sedov blast the relative energy error is :math:`3 \times 10^{-3}` with the
  hottest-neighbour sound speed at 0.2 and 0.16 to 0.23 without it. The
  time-step limiter does not repair this, since the acceleration criterion
  already keeps the hot particles' neighbours within its reach. The
  hottest-neighbour term is therefore kept.
+ The step count of this scheme at :math:`C_{\rm CFL} = 0.2` equals that of
  SPHENIX at 0.1, the two schemes' recommended values, at equal accuracy:
  the factor two between the conventions is not a cost.

Kernel and neighbour number
---------------------------

The paper uses the Wendland C6 kernel with 300 neighbours
(:math:`\eta = 1.70` in SWIFT's convention). The recommended configuration
in SWIFT is the quintic spline with ``SPH:resolution_eta: 1.5`` (150
neighbours in 3D, 33 in 2D). In the 2026-10 kernel comparison on the test
suite at :math:`\eta = 1.5` the quintic gave the most accurate densities and
shock profiles at the lowest cost; the cubic (86 neighbours) and quartic
(116) splines pair at this neighbour number, and the Wendland kernels carry
a positive density bias of 0.3-2 per cent at these neighbour numbers (the
self-contribution bias of Dehnen & Aly 2012, not corrected in SWIFT).
Wendland C6 (206 neighbours) is the alternative when the smoothest velocity
field in a subsonic flow matters more than cost. The scheme does not change
SWIFT's compile-time default (cubic spline), so ``--with-kernel=quintic-spline``
must be given explicitly, and the generic example parameter files must be run
with ``-P SPH:resolution_eta:1.5 -P SPH:CFL_condition:0.2``.

Departures from Rosswog (2020)
------------------------------

The choices below change results with respect to the paper and are
deliberate; they are listed so that comparisons with the paper can account
for them.

+ **Viscosity in SWIFT's** :math:`h`. The paper's kernels have a support of
  :math:`2h`, SWIFT's of :math:`\gamma_K h`, and :math:`\mu_a` is linear in
  :math:`h`. The limiter's lengths are converted; :math:`\mu_a` is not, so at
  the same kernel support the code's viscosity is the paper's with
  :math:`\alpha` and :math:`\beta` scaled by :math:`2/\gamma_K` and
  :math:`(2/\gamma_K)^2`: with the quintic spline, :math:`\alpha = 1`,
  :math:`\beta = 2` act as the paper's 0.91 and 1.66, with the paper's own
  Wendland C6 as 0.82 and 1.33. The values were tuned on the test suite as
  they stand and are kept.
+ **CFL factor of two**, see the time step above.
+ **Linear reconstruction** instead of the paper's quadratic one, by
  design. The paper reports a moderate effect (sharper interfaces, faster
  Kelvin-Helmholtz growth).
+ **Zeroth-order force error.** The matrix correction makes difference
  gradients linearly exact, but the pair operator :math:`\mathbf{G}_a` does
  not sum to zero over the neighbours, so on disordered particles the
  residual pressure force at constant pressure is about 1.5 times that of
  plain SPH (CRK-SPH's symmetrised form has 0.45 times). This is inherent to
  the formulation, not to the implementation.
+ **Conduction signal speed.** The paper's eq. 25 gives a pressure-based
  signal speed without gravity and the velocity-based eq. 26 with gravity;
  the code always uses the velocity-based form (plus the pressure-residual
  term above), chosen at compile time in ``hydro_parameters.h``.
+ **Dissipation switches not in the paper**: the time-dependent viscosity
  coefficient in place of the constant :math:`\alpha = 1` (restored with
  ``SPH:viscosity_alpha_min: 1``), the Cullen & Dehnen shear
  limiter on :math:`Q`, with its Mach threshold and dead zone, the
  convergence-scaled cut-off of the reconstruction, the pressure-residual
  conduction term and the per-particle conduction coefficient with its
  viscosity cap (the paper's conduction is on in shocks, the code's is
  off), the gates that apply :math:`Q` only to approaching pairs with
  positive heating, the neighbour-maximum sound speed in the time step, and
  SWIFT's leapfrog with individual time steps instead of the paper's global
  TVD Runge-Kutta step.

Parameters
----------

Most of the scheme's constants are fixed at compile time in
``src/hydro/MAGMA/hydro_parameters.h`` (changing one means editing that
file and recompiling) at the values settled by the
2026-10 tuning study on the test suite (:math:`\alpha_{\max} = 1`,
:math:`\beta = 2`, :math:`\epsilon = 0.1`, the cut-off width, the shear
limiter and its Mach threshold 0.5 and dead zone, the angle limit 0.5 rad,
:math:`\alpha_u = 0.05`, :math:`\alpha_P = 1`, :math:`f_0 = 0.05`). The
run-time parameters and their defaults are:

.. code-block:: yaml

    SPH:
      resolution_eta:                 1.5    # Recommended (compulsory parameter, no default)
      CFL_condition:                  0.2    # Recommended (compulsory parameter, no default)

      limiter_eta_crit:              -1.     # Cut-off distance of the reconstruction in units of h; < 0: 1 / resolution_eta (the paper's value)
      gradient_max_condition_number: 60.     # Clipping of the inverse C-matrix in degenerate directions
      viscosity_limiter_floor:        0.     # xi = max(xi, floor): fraction of the viscosity kept in pure shear
      viscosity_alpha_min:            0.5    # Floor of the per-particle viscosity alpha (its maximum is 1); 1 = the paper's constant viscosity
      viscosity_length:               0.05   # Cullen & Dehnen's l: its decay time is H / (2 l c_s)
      diffusion_beta:                 1.0    # Growth rate of the per-particle conduction coefficient (SPHENIX's beta_D)
      diffusion_alpha_min:            1.0    # Floor (and initial value) of the per-particle conduction coefficient; 0 = SPHENIX's behaviour

The generic keys ``viscosity_alpha``, ``viscosity_alpha_max``,
``viscosity_beta``, ``diffusion_alpha`` and ``diffusion_alpha_max`` are not
read by this scheme. The parameters are validated at start-up
(``viscosity_init()`` and ``diffusion_init()``): the run aborts for a
condition number below 1, a floor outside [0, 1], ``viscosity_alpha_min``
above 1 or ``diffusion_alpha_min`` above 1, or a non-positive
``viscosity_length``.

Snapshot output
---------------

Besides the standard fields, the scheme writes per particle:
``GradientsXVelocities``, ``GradientsYVelocities``, ``GradientsZVelocities``
and ``GradientsInternalEnergies`` (the matrix gradients, comoving),
``FallbackFlags`` (bit 0: C-matrix not invertible, bit 1: :math:`h` near
``h_max``, both meaning standard SPH gradients were used; bit 2: the inverse
was regularised but the gradient functions were still used),
``GradientFallbackPairs`` (pairs that used kernel gradients in the last force
loop), ``ConvergenceFraction``, ``ViscosityLimiter`` (the Cullen & Dehnen
:math:`\xi_a`), ``DiffusionParameters`` (:math:`\alpha^D_a`),
``ViscosityParameters`` (the per-particle :math:`\alpha_a`; the names are
SPHENIX's), ``VelocityDivergences`` (the trace of the matrix velocity
gradient, comoving, without the Hubble flow) and
``LaplacianInternalEnergies`` (the SPH Laplacian of :math:`u` of the last
gradient loop). A healthy run has ``FallbackFlags`` equal to 0 for essentially all
particles and a mean ``GradientFallbackPairs`` well below 1; the
scheme-comparison suite reports both fractions in its summary table. The
``HydroScheme`` group records every constant and parameter of the
dissipation terms (``Alpha viscosity``, ``Beta viscosity``,
``Epsilon viscosity``, ``Limiter eta_crit``, ``Limiter width``,
``Viscosity limiter floor``,
``Viscosity limiter Mach``, ``Viscosity boost dead zone (noise units)``,
``Alpha viscosity (min)``, ``Viscosity decay length``,
``Gradient max condition number``, ``Gradient cos angle limit``,
``Alpha diffusion``, ``Alpha diffusion (pressure residual)``,
``Pressure residual threshold``, ``Beta diffusion``,
``Alpha diffusion (min)``).

Feedback and sub-grid hooks
---------------------------

The scheme implements the pressure floor and the entropy floor, the
minimal internal energy, the time-step limiter, the velocity-kick signal
velocity correction of kinetic feedback (the Courant signal speed is raised
by :math:`\beta\,\Delta v`, as in the other schemes) and the diffusive
feedback reset (per-particle viscosity coefficient to its maximum and
conduction coefficient to zero until the particle's next force preparation;
the raised coefficient, through the conduction cap, keeps the conduction of
the particle and its neighbours off until it has decayed, about ten
sound-crossing times of the kernel support, whatever the viscosity floor).
The velocity
divergence used by the GEAR star formation and sink models is the trace of
the matrix velocity gradient recorded at the last force preparation (zero
for particles that fell back to standard SPH). Cosmological runs are
supported: the Hubble flow enters the approach test, :math:`\mu_a`, the
viscous heating, the conduction signal speed, the velocity limiter and the
divergence seen by the shear limiter, the shock-presence sum and the
:math:`\alpha_a` trigger, and the adiabatic :math:`P\,dV` term uses the
peculiar velocity only.

Known limitations
-----------------

+ Adaptive softening is not implemented (configuring with
  ``--with-adaptive-softening`` fails to compile).
+ The reconstruction is linear; together with the gradient noise this
  leaves a few per cent of the velocity jump to the viscosity in smooth
  flow, so a linear sound wave is damped by 2 per cent per period at the
  suite's resolutions with the default floor of the viscosity coefficient
  (see the ``SoundWave_3D`` and ``WaveDamping`` tests of the
  scheme-comparison suite; SPHENIX damps it by 0.5 per cent).
+ The scheme is not hand-vectorised (``--enable-hand-vec`` applies to
  GADGET-2 only).
+ Angular momentum is conserved only approximately (as in the paper).
+ Particles at ``h_max`` and particles whose correction matrix cannot be
  inverted use plain SPH without reconstruction.
+ Three loops over the neighbours, a 3x3 inversion per particle and the
  reconstruction in every force pair make a step about 2.5 to 3 times the
  cost of a SPHENIX step at the same neighbour number, and the recommended
  neighbour number is larger than SPHENIX's.
+ MHD, radiative transfer and the GIZMO-style Riemann options have not been
  tested with this scheme.

Examples
--------

The scheme is exercised by ``examples/HydroTests/SchemeComparison``
(``run_suite.sh``, which builds MI1, MI2 and MI3 beside SPHENIX from the same
tree and runs 27 standard tests with identical initial conditions) and by
the individual hydro tests, in particular those added for the scheme's
dissipation switches: ``KeplerianRing_CullenDehnen`` (shear limiter and Mach
threshold), ``KeplerianRing`` in 2D (the cold Mach-700 disc the dead zone
was set on), ``ShearingShock_3D`` (a shock seen through strong shear),
``ThermalDiffusion_2D`` (pressure-residual conduction at a contact),
``SoundWave_3D`` and ``TaylorGreenVortex_2D`` (numerical dissipation in
smooth flow), ``StrongBlastWave_3D`` (section 3.4.3 of the paper),
``YeeVortex_2D`` (order of convergence) and ``DrivenTurbulence_3D``. In each
case run the example with ``-P SPH:resolution_eta:1.5 -P SPH:CFL_condition:0.2``
on top of its parameter file, or edit the file.

What to expect at the recommended configuration: Sod and Sedov profiles at
least as sharp as SPHENIX's with a smaller velocity error; the square test
keeps its corners (MI2/MI3) or becomes a diamond (MI1); the Gresho vortex,
the Cullen & Dehnen ring and the cold 2D disc survive with the shear
limiter and spread within a few orbits without it; the Yee vortex converges
at first order in :math:`L_1(\rho)` with a slope of about :math:`-0.8`.

References
----------

+ S. Rosswog, "The Lagrangian hydrodynamics code MAGMA2", MNRAS 498, 4230
  (2020), `arXiv:1911.13093 <https://arxiv.org/abs/1911.13093>`_: the
  scheme, its MI1 and MI2 formulations, the slope-limited reconstruction and
  the dissipation terms.
+ S. Rosswog, "Boosting the accuracy of SPH techniques: Newtonian and
  special-relativistic tests", MNRAS 448, 3628 (2015): the matrix-inversion
  gradients.
+ L. Cullen & W. Dehnen, "Inviscid SPH", MNRAS 408, 669 (2010),
  `arXiv:1006.1524 <https://arxiv.org/abs/1006.1524>`_: the shear limiter
  and the per-particle viscosity coefficient.
+ N. Frontiere, C. D. Raskin & J. M. Owen, "CRKSPH", J. Comp. Phys. 332, 160
  (2017), `arXiv:1605.00725 <https://arxiv.org/abs/1605.00725>`_: the
  van Leer limiter of the reconstruction.
+ D. J. Price, "Modelling discontinuities and Kelvin-Helmholtz
  instabilities in SPH", J. Comp. Phys. 227, 10040 (2008): the
  pressure-jump conduction signal speed.
+ J. Borrow, M. Schaller, R. G. Bower & J. Schaye, "SPHENIX", MNRAS 511,
  2367 (2022): the per-particle conduction coefficient and its viscosity
  cap, the feedback reset.
+ W. Dehnen & H. Aly, "Improving convergence in SPH simulations without
  pairing instability", MNRAS 425, 1068 (2012): kernel choice and
  neighbour numbers.
