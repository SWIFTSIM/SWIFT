.. SIDM
   Katy Proctor, August 2025

.. _SIDM:

Self-interacting dark matter (SIDM) implementation
==================================================

.. warning::
    The SIDM scheme is still in development and is not useable
    at this moment. This page is currently a placeholder to document new
    features and requirements as the code grows.

SIDM is implemented as a separate particle species using ``PartType7`` (i.e.
the code does not consider ``PartType1`` particles as being of SIDM nature).

Each SIDM particle has its own smoothing length, found with a neighbour loop
like the one used for the gas. Pairs of neighbouring particles then scatter
with a probability set by the overlap of their two kernels, following the
kernel-overlap approach of `Fischer et al. (2026)
<https://scixplorer.org/abs/2026OJAp....971159F/abstract>`_.


Configuring and running
-----------------------

SIDM is switched on at configure time and needs three choices:

.. tabularcolumns:: |p{0.33\linewidth}|p{0.38\linewidth}|p{0.19\linewidth}|

.. list-table::
   :header-rows: 1
   :widths: 30 45 25

   * - Configure option
     - Meaning
     - Default
   * - ``--with-sidm``
     - SIDM model: ``none`` or ``Basic``.
     - ``none``
   * - ``--with-sidm-kernel``
     - Kernel for the SIDM neighbour finding loop: ``cubic-spline``, ``quartic-spline``,
       ``quintic-spline``, ``wendland-C2``, ``wendland-C4`` or
       ``wendland-C6``.
     - Same as ``--with-kernel`` (``cubic-spline`` by default)
   * - ``--with-sidm-overlap-kernel``
     - Kernel used for the pair overlap in the scattering loop:
       ``cubic-spline`` or ``tophat``.
     - ``cubic-spline``

At run time, SIDM is enabled with the ``--sidm`` flag::

  ./swift --sidm --self-gravity --threads=8 params.yml


Parameters
----------

The parameters are read from the ``BasicSIDM`` section of the parameter file:

.. code:: YAML

   BasicSIDM:
     resolution_eta:         1.2348  # Smoothing length in units of the mean inter-particle separation (48 neighbours with the cubic spline).
     constant_cross_section: 1       # Use a constant cross-section per unit mass.
     sigma_over_m_cgs:       1.      # Cross-section per unit mass [cm^2/g].
     h_tolerance:            1e-4    # (Optional) Relative accuracy of the smoothing-length iteration.
     h_max:                  10.     # (Optional) Maximal smoothing length in internal units.
     h_min_ratio:            0.      # (Optional) Minimal smoothing length in units of the DM softening (divided by the kernel gamma).
     max_ghost_iterations:   30      # (Optional) Maximal number of smoothing-length iterations.
     max_volume_change:      1.4     # (Optional) Maximal change of the kernel volume over one time-step.
     kappa_timestep:         0.01    # (Optional) Time-step criterion: dt < kappa_timestep / Gamma.

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Parameter
     - Default
     - Notes
   * - ``resolution_eta``
     - required
     - :math:`\eta` in :math:`h = \eta\,(m/\rho)^{1/3}`.
   * - ``constant_cross_section``
     - required
     - Only ``1`` (constant :math:`\sigma/m`) is implemented.
   * - ``sigma_over_m_cgs``
     - required
     - Read when ``constant_cross_section`` is ``1``.
   * - ``h_tolerance``
     - ``1e-4``
     -
   * - ``h_max``
     - no limit
     -
   * - ``h_min_ratio``
     - ``0``
     -
   * - ``max_ghost_iterations``
     - ``30``
     -
   * - ``max_volume_change``
     - hydro value
     - Uses the hydro scheme's value if not given.
   * - ``kappa_timestep``
     - ``0.01``
     -


Neighbour finding
-----------------

The smoothing length :math:`h_i` of each SIDM particle is found with the same
iterative scheme as the gas (`Price 2012
<https://ui.adsabs.harvard.edu/abs/2012JCoPh.231..759P>`_). The SIDM ghost
adjusts :math:`h_i` until the kernel density, which includes the particle's own
contribution, satisfies

.. math::
   \rho_i = \sum_j m_j W(|\mathbf{x}_i - \mathbf{x}_j|, h_i), \qquad
   h_i = \eta \left(\frac{m_i}{\rho_i}\right)^{1/3},

to a relative accuracy ``h_tolerance``. Neighbours are searched for within
:math:`H_i = \gamma h_i`, where :math:`\gamma` is the support of the chosen
density kernel in units of :math:`h` (table 1 of `Dehnen & Aly 2012
<https://ui.adsabs.harvard.edu/abs/2012MNRAS.425.1068D>`_). The corresponding
neighbour number is

.. math::
   N_{\rm ngb} = \frac{4\pi}{3} (\gamma \eta)^3 .

The same :math:`\eta` therefore gives different neighbour numbers for
different kernels:

.. list-table::
   :header-rows: 1

   * - ``--with-sidm-kernel``
     - :math:`\gamma`
     - :math:`N_{\rm ngb}` at :math:`\eta = 1.2348`
     - :math:`\eta` for 48 neighbours
   * - ``cubic-spline``
     - 1.825742
     - 48.0
     - 1.2348
   * - ``quartic-spline``
     - 2.018932
     - 64.9
     - 1.1167
   * - ``quintic-spline``
     - 2.195775
     - 83.5
     - 1.0267
   * - ``wendland-C2``
     - 1.936492
     - 57.3
     - 1.1642
   * - ``wendland-C4``
     - 2.207940
     - 84.9
     - 1.0211
   * - ``wendland-C6``
     - 2.449490
     - 115.9
     - 0.9204


Scattering
----------

In the force loop, each pair of neighbours :math:`i, j` (separation
:math:`r_{ij} < \gamma\,\max(h_i, h_j)`) contributes a pair scattering rate

.. math::
   \Gamma_{ij} = \tfrac{1}{2}(m_i + m_j)\,\frac{\sigma}{m}\,|\mathbf{v}_i - \mathbf{v}_j|\,\Lambda_{ij},

where :math:`\Lambda_{ij}` is the kernel overlap described below. Velocities
and :math:`\Lambda_{ij}` are converted to physical units, so in cosmological
runs :math:`\Lambda_{ij}` carries a factor :math:`a^{-3}`. Each particle's
total rate is :math:`\Gamma_i = \sum_{j \neq i} \Gamma_{ij}` (written to the
snapshots as ``Rates``).

The pair scatters with probability

.. math::
   P_{ij} = \Gamma_{ij}\,\min(\Delta t_i, \Delta t_j).

The shorter time-step is used because the earlier part of the longer time-step has
already been covered by the non-symmetric interactions computed while only the
particle on the shorter time-step was active. Scattering is elastic and isotropic
in the centre-of-mass frame.

The time-step of each particle is limited to
:math:`\Delta t_i < \kappa / \Gamma_i`, with :math:`\kappa` =
``kappa_timestep``, so that the scattering probability per step stays small.


Overlap kernels
---------------

The overlap of the two particles' kernels is

.. math::
   \Lambda_{ij} = \int W(|\mathbf{x} - \mathbf{x}_i|, h_i)\,
                       W(|\mathbf{x} - \mathbf{x}_j|, h_j)\,{\rm d}^3x .

Because each kernel integrates to one, so does :math:`\Lambda_{ij}` over all
separations. It depends only on :math:`r = r_{ij} / h_{\max}` and
:math:`s = h_{\min} / h_{\max}`, up to a factor :math:`h_{\max}^{-3}`.

**Truncation.** The full overlap extends to :math:`r_{ij}` beyond the
neighbour search radius :math:`\gamma h_{\max}`, so those pairs are never
seen by the force loop. Both overlap kernels are therefore truncated at
:math:`\gamma h_{\max}` and divided by

.. math::
   N(s) = \int_0^{\gamma} 4\pi r^2\, \Lambda(r, s)\, {\rm d}r
   \qquad (h_{\max} = 1),

so that the truncated overlap still integrates to one. Here :math:`\gamma` is
the support of the *density* kernel.

Top-hat
~~~~~~~

Each particle is a uniform sphere of radius :math:`h` (not :math:`\gamma h`)
with :math:`W = 1/V`, so :math:`\Lambda_{ij} = V_{\rm lens}/(V_i V_j)`, where
the intersection volume of the two spheres is (`Weisstein, "Sphere-Sphere
Intersection" <https://mathworld.wolfram.com/Sphere-SphereIntersection.html>`_)

.. math::
   V_{\rm lens} = \frac{\pi (h_i + h_j - r)^2
                  \left[r^2 + 2r(h_i + h_j) - 3(h_i - h_j)^2\right]}{12\,r}

for :math:`|h_i - h_j| < r < h_i + h_j`, and
:math:`V_{\rm lens} = \tfrac{4\pi}{3}h_{\min}^3` when the smaller sphere lies
inside the larger one. The normalisation has a closed form,
:math:`N = 1 - M` with

.. math::
   M = \frac{3}{16 s^3}\,T^3 \left[4as - (a^2 + 3s)\,T + aT^2 - \frac{T^3}{6}\right],
   \qquad a = 1 + s, \quad T = a - \gamma,

and :math:`M = 0` when :math:`T \le 0` (the lens never reaches the search
radius).

Cubic spline
~~~~~~~~~~~~

For two cubic splines the overlap is computed by integrating over spherical
shells of radius :math:`\rho = s t` around the particle with the smaller
kernel. The average of the other kernel over such a shell is
:math:`[G(r + \rho) - G(|r - \rho|)] / (2 r \rho)`, with
:math:`G(b) = \int_0^b u\,W(u)\,{\rm d}u` known in closed form. This gives
the one-dimensional integral

.. math::
   \Lambda(r, s) = \frac{2\pi}{r s} \int_0^{\gamma_{\rm cs}} t\,W(t)
                   \left[G(r + s t) - G(|r - s t|)\right]{\rm d}t ,

with :math:`\gamma_{\rm cs} = 1.825742` and :math:`h_{\max} = 1`. The limits
are :math:`\Lambda = 4\pi\int t^2 W(t) W(st)\,{\rm d}t` for :math:`r \to 0` and
:math:`\Lambda = W(r)` for :math:`s \to 0`.

:math:`\Lambda` and :math:`1/N(s)` are tabulated once at start-up on a grid of 128 values of :math:`s` and 256 values of
:math:`x = r_{ij} / [\gamma_{\rm cs}(h_i + h_j)]`, using Simpson's rule with
512 intervals, and read by bilinear interpolation. The integrals agree with
direct 3D quadrature to :math:`\sim 10^{-10}`.
The code is in ``src/kernel_sidm.h`` and ``src/kernel_sidm.c``.