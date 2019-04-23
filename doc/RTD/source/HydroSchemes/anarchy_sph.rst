.. ANARCHY-SPH
   Josh Borrow 5th April 2018

ANARCHY-PU SPH
==============

This scheme is similar to the one used in the EAGLE code. This scheme
includes:

+ Durier & Dalla Vecchia (2012) time-step limiter
+ Pressure-Energy SPH
+ Thermal diffusion following Price (2012)
+ A simplified version of the 'Inviscid SPH' artificial viscosity
  (Cullen & Denhen 2010).

More information will be made available in a forthcoming publication.

The scheme as-implemented in SWIFT is slightly different to the one
implemented in the original EAGLE code:

+ Pressure-Energy SPH is used instead of Pressure-Entropy SPH
+ Artificial viscosity coefficients have changed -- from minimal
  value of 0.1 to 0.0, and from length of 0.1 to 0.25. This
  is based on performance of hydrodynamics tests in SWIFT and may
  be to do with our choice of smoothing length definition.
+ Recommended kernel changed from Wendland-C2 (with 100 Ngb) to
  Quintic Spline (with ~82 Ngb).


.. code-block:: bash
   
   ./configure --with-hydro=anarchy-pu --with-kernel=quintic-spline --disable-hand-vec

