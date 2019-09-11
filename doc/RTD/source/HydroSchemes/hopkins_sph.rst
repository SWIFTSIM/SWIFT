.. 'Hopkins'-SPH
   Josh Borrow 5th April 2018

Pressure-Entropy SPH
====================

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Contents:

A Pressure-Entropy SPH scheme is available in SWIFT, inspired by Hopkins 2013.
This includes a fixed Monaghan AV scheme and a Balsara switch.


.. code-block:: bash
   
   ./configure --with-hydro="pressure-entropy"


The parameters available for this scheme, and their defaults, are:

.. code-block:: yaml

    SPH:
        viscosity_alpha: 0.8  # Fixed value for the alpha viscosity



Pressure-Energy SPH
===================

Pressure-energy SPH is now implemented in SWIFT, and like the pressure-entropy
scheme it includes a Monaghan AV scheme and a Balsara switch.


.. code-block:: bash
   
   ./configure --with-hydro="pressure-energy"


The parameters available for this scheme, and their defaults, are:

.. code-block:: yaml

    SPH:
        viscosity_alpha: 0.8  # Fixed value for the alpha viscosity


There is a variant of this implementation that includes a Morris & Monaghan
(1997) variable artificial viscosity that aims to reduce disappation away
from strong shocks. This implementation also includes a Balsara switch.
To use this scheme, you should use:

.. code-block:: bash
   
   ./configure --with-hydro="pressure-energy-monaghan"


The parameters available for this scheme, and their defaults, are:

.. code-block:: yaml

    SPH:
        viscosity_alpha: 0.8  # Initial value for the alpha viscosity
        viscosity_length: 0.25  # Viscosity decay length (in terms of sound-crossing time)
        # These are enforced each time-step
        viscosity_alpha_max: 2.0  # Maximal allowed value for the viscosity alpha
        viscosity_alpha_min: 0.1  # Minimal allowed value for the viscosity alpha


There is also a compile-time parameter, ``viscosity_beta`` that we set to
3.0. During feedback events, the viscosity is set to the compile-time
``hydro_props_default_viscosity_alpha_feedback_reset = 2.0``. These can be
changed in ``src/hydro/AnarchyPU/hydro_parameters.h``.
