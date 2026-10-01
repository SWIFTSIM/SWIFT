.. SPH kernels
   Matthieu Schaller, October 2026

.. _sph_kernels:

SPH kernels
===========

Six kernels are available via the ``--with-kernel`` configuration option:
``cubic-spline`` (default), ``quartic-spline``, ``quintic-spline``,
``wendland-C2``, ``wendland-C4`` and ``wendland-C6``. The Wendland C4 and C6
kernels are not defined in 1D.

All constants follow table 1 of `Dehnen & Aly (2012)
<https://ui.adsabs.harvard.edu/abs/2012MNRAS.425.1068D>`_. The kernel is

.. math::

   W(r, h) = \frac{C}{H^d} \, f\!\left(\frac{r}{H}\right), \qquad H = \gamma h,

with :math:`d` the dimension, :math:`H` the compact support and
:math:`\gamma` chosen such that the number of neighbours implied by
``SPH:resolution_eta`` (see :ref:`Parameters_SPH`) is the same for all
kernels. In the code ``kernel_gamma`` is :math:`\gamma` and the functions
``kernel_eval()`` and ``kernel_deval()`` take :math:`u = r/h` as argument.

Evaluation
----------

The kernels are stored in ``src/kernel_hydro.h`` as piecewise polynomials
in :math:`t = u - u_0` on uniform sub-intervals of :math:`[0, \gamma)`,
with the normalisation folded into the coefficients and the origin
:math:`u_0` of each sub-interval placed at its right end. The last origin is
the edge of the support, so :math:`W(\gamma h) = 0` exactly and the relative
accuracy of :math:`W` and :math:`\partial W / \partial r` stays at a few
:math:`10^{-7}` all the way to the edge. (A polynomial in :math:`r/H`
expanded about 0 behaves as :math:`(1 - r/H)^n` near the edge and loses all
relative accuracy there by cancellation.) The constants are computed in
double precision from exact expressions and rounded once. The
double-precision version ``kernel_eval_double()`` uses the single-precision
:math:`\gamma` promoted to double, so that both versions have the same
support.

The hand-vectorised functions used by the GADGET-2 scheme use the same
tables and the same arithmetic as the scalar ones.

The accuracy is checked against the analytic expressions by
``tests/testKernelAccuracy`` for every kernel and dimension.
