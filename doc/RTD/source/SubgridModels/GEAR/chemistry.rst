.. GEAR sub-grid model chemistry
   Darwin Roduit, 30th March 2025

.. gear_chemistry:

.. _gear_chemistry:

Chemistry
=========

Feedback mechanisms such as supernova feedback transfer metals to the nearby gas particles. There is no mass exchange (mass advection) in SPH methods, as in grid-based or moving-mesh hydrodynamics solvers. As a result, the metals received by the gas particles are locked in these particles. Grid or moving-mesh codes have mass advection and often implement passive scalar advection of metals to model metal mixing. GEAR implements different methods to model metal mixing.

.. _gear_smoothed_metallicity:

Smoothed metallicity
--------------------

The smoothed metallicity scheme consists in using the SPH to smooth the metallicity of each particle over the neighbors. It is worth to point the fact that we are *not exchanging* any metals but only smoothing it. The parameter ``GEARChemistry:initial_metallicity`` set the (non smoothed) initial mass fraction of each element for all the particles and ``GEARChemistry:scale_initial_metallicity`` use the feedback table to scale the initial metallicity of each element according the Sun's composition. If ``GEARChemistry:initial_metallicity`` is negative, then the metallicities are read from the initial conditions.

For this chemistry scheme the parameters are:

.. code:: YAML

   GEARChemistry:
    initial_metallicity: 1         # Initial metallicity of the gas (mass fraction)
    scale_initial_metallicity: 1   # Should we scale the initial metallicity with the solar one?

