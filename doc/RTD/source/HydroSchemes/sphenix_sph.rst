.. SPHENIX SPH
   Josh Borrow 8th January 2020

SPHENIX
=======

This is the new scheme that will be used in EAGLE-XL. This scheme includes:

+ Durier & Dalla Vecchia (2012) time-step limiter
+ Density-Energy SPH
+ Thermal diffusion following Price (2012)
+ A simplified version of the 'Inviscid SPH' artificial viscosity
  (Cullen & Dehnen 2010), with a Balsara switch
+ A diffusion limiter, used to prevent energy leakage out of EAGLE
  supernovae (Borrow+ 2020).


The simplified version of the 'Inviscid SPH' artificial viscosity calculates
the time differential of the velocity divergence explicitly, using the value
from the previous step. We also use the Balsara switch instead of the improved
neighbour-based limiter from Cullen & Dehnen 2010, to avoid matrix
calculations.


To configure with this scheme, use

.. code-block:: bash
   
   ./configure --with-hydro=sphenix --with-kernel=quintic-spline --disable-hand-vec


The diffusion limiter is implemented to ensure that the diffusion is turned
off in very viscous flows and works as follows:

.. code-block:: C

    float new_diffusion_alpha = old_diffusion_alpha;

    const float viscous_diffusion_limit =
      diffusion_alpha_max *
      (1.f - maximum_alpha_visc_over_ngb / viscosity_alpha_max);

    new_diffusion_alpha = min(new_diffusion_alpha, viscous_diffusion_limit);


The parameters available for this scheme, and their defaults, are:

.. code-block:: yaml

    SPH:
        viscosity_alpha: 0.1  # Initial value for the alpha viscosity
        viscosity_length: 0.25  # Viscosity decay length (in terms of sound-crossing time)
        # These are enforced each time-step
        viscosity_alpha_max: 2.0  # Maximal allowed value for the viscosity alpha
        viscosity_alpha_min: 0.0  # Minimal allowed value for the viscosity alpha

        diffusion_alpha: 0.0  # Initial value for the diffusion alpha
        diffusion_beta: 0.25  # Timescale to raise the diffusion coefficient over
                              # (decay is on the sound-crossing time)
        # These are enforced each time-step
        diffusion_alpha_max: 1.0
        diffusion_alpha_min: 0.0


There is also a compile-time parameter, ``viscosity_beta`` that we set to
3.0. During feedback events, the viscosity is set to the compile-time
``hydro_props_default_viscosity_alpha_feedback_reset = 2.0`` and the
diffusion is set to ``hydro_props_default_diffusion_alpha_feedback_reset =
0.0``. These can be changed in ``src/hydro/SPHENIX/hydro_parameters.h``.
