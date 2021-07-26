.. GASOLINE-SPH
   Josh Borrow 27 March 2021

Gasoline-2 SPH
==============

This scheme very closely follows the paper by Wadsley et al. (2017).
It includes a "Geometric Density Average Force" equation of motion,
explicit corrections for gradients, a Cullen & Dehnen (2010) style
artificial viscosity and switch, and an artificial conduction
scheme based on the local shear tensor.

The equation of motion:

.. math::

   \frac{\mathrm{d} \mathbf{v}_{i}}{\mathrm{d} t}=-\sum_{j} m_{j}\left(\frac{P_{i}+P_{j}}{\rho_{i} \rho_{j}}+\Pi_{i j}\right) \nabla_{i} \bar{W}_{i j}

.. math::
   
   \frac{\mathrm{d} u_{i}}{\mathrm{d} t}=\sum_{j} m_{j}\left(\frac{P_{i}}{\rho_{i} \rho_{j}}+\frac{1}{2} \Pi_{i j}\right) \mathbf{v}_{i j} \cdot \nabla_{i} \bar{W}_{i j}

with

.. math::

   \nabla_{i} \bar{W}_{i j}=\frac{1}{2} f_{i} \nabla_{i} W\left(r_{i j}, h_{i}\right)+\frac{1}{2} f_{j} \nabla_{j} W\left(r_{i j}, h_{j}\right)

.. math::
   
   f_{i}=\frac{\sum_{j} \frac{m_{j}}{\rho_{i}} \mathbf{r}_{i j}^{2} W^{\prime}\left(\frac{r_{i j}}{h_{i}}\right)}{\sum_{j} \frac{m_{j}}{\rho_{j}} \mathbf{r}_{i j}^{2} W^{\prime}\left(\frac{r_{i j}}{h_{i}}\right)}


For more information, please see the Wadsley et al. (2017) paper. Note that we do not
use the exact same time-stepping condition, but instead use the same as the other
schemes in SWIFT.

To configure with this scheme, use

.. code-block:: bash
   
   ./configure --with-hydro=gasoline --with-kernel=wendland-C4 --disable-hand-vec


The parameters available for this scheme, and their defaults, are:

.. code-block:: yaml

    SPH:
        viscosity_alpha: 0.1  # Initial value for the alpha viscosity
        viscosity_length: 0.2  # Viscosity decay length (in terms of sound-crossing time)
        # These are enforced each time-step
        viscosity_alpha_max: 2.0  # Maximal allowed value for the viscosity alpha
        viscosity_alpha_min: 0.0  # Minimal allowed value for the viscosity alpha

        diffusion_coefficient: 0.03 # Controls the diffusion rate


There is also a compile-time parameter, ``viscosity_beta`` that we set to
3.0. During feedback events, the viscosity is set to the compile-time
``hydro_props_default_viscosity_alpha_feedback_reset = 2.0`` and the
diffusion is set to ``hydro_props_default_diffusion_rate_feedback_reset =
0.0``. These can be changed in ``src/hydro/Gasoline/hydro_parameters.h``.

