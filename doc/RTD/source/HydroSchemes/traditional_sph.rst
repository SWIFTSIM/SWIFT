.. Traditional SPH (GADGET-2)
   Josh Borrow 4th April 2018

Traditional (Density-Entropy) SPH
=================================

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Contents:

Traditional, GADGET-2-like, Density-Entropy SPH is available in SWIFT with
a Monaghan artificial viscosity scheme and Balsara switch.

To use this hydro scheme, you need no extra configuration options -- it is the
default!

As it uses a very simple, fixed artificial viscosity, only the
``SPH:viscosity_alpha`` parameter has any effect for this scheme. This will
change the strength of the artificial viscosity throughout the simulation,
and has a default of 0.8.

