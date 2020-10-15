.. PHANTOM SPH
   Josh Borrow 13th October 2020

Phantom
=======

This scheme is a reference implementation similar to the one presented in
Price (2018), the PHANTOM paper (not including MHD). It uses:

+ A simplified Cullen & Dehnen AV limiter (note that this is different to 
  PHANTOM as we do not explicitly include the matrix calculation).
+ A fixed alpha artificial conduction scheme used for hydro-only problems
  as presented in the PHANTOM paper (i.e. we use the 'hydro-only' conduction
  velocity, rather than the one used for gravitational problems).
+ Base Density-Energy SPH

The simplified version of the 'Inviscid SPH' artificial viscosity calculates
the time differential of the velocity divergence explicitly, using the value
from the previous step. We also use the Balsara switch instead of the improved
neighbour-based limiter from Cullen & Dehnen 2010, to avoid matrix
calculations. We also use a different value for the 'h-factors' due to SWIFT
using neighbour finding based on particle number density, rather than local
mass density.


To configure with this scheme, use

.. code-block:: bash
   
   ./configure --with-hydro=phantom --disable-hand-vec


The parameters available for this scheme, and their defaults, are:

.. code-block:: yaml

    SPH:
        viscosity_alpha: 0.1  # Initial value for the alpha viscosity
        viscosity_length: 0.25  # Viscosity decay length (in terms of sound-crossing time)
        # These are enforced each time-step
        viscosity_alpha_max: 2.0  # Maximal allowed value for the viscosity alpha
        viscosity_alpha_min: 0.0  # Minimal allowed value for the viscosity alpha

        diffusion_alpha: 1.0  # Fixed value for the diffusion alpha


There is also a compile-time parameter, ``viscosity_beta`` that we set to
3.0. During feedback events, the viscosity is set to the compile-time
``hydro_props_default_viscosity_alpha_feedback_reset = 2.0`` and the
diffusion is set to ``hydro_props_default_diffusion_alpha_feedback_reset =
0.0``. These can be changed in ``src/hydro/Phantom/hydro_parameters.h``.

