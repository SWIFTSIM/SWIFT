.. AGORA sub-grid model
   Yves Revaz 04.2023

.. _agora:

AGORA model
===========

AGORA's model corresponds to the subgrid model adopted by the 
`AGORA Project <https://sites.google.com/site/santacruzcomparisonproject/>`_ (High-resolution Galaxy Simulations Comparison Project).
Details of the model is described in 
`Kim \& al. 2014 <https://ui.adsabs.harvard.edu/link_gateway/2014ApJS..210...14K/PUB_PDF>`_
and
`Kim \& al. 2016 <https://ui.adsabs.harvard.edu/link_gateway/2016ApJ...833..202K/PUB_PDF>`_.



This model can be selected with the configuration option ``--with-subgrid=AGORA`` and run with the option ``--agora``. 
A few examples exist and can be found in ``examples/AGORA``. 

.. _agora_grackle_cooling:

Gas cooling/heating: Grackle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The AGORA model uses the Grackle cooling library which is also used by the :ref:`GEAR model <gear_grackle_cooling>`.




.. _agora_star_formation:

Star formation
~~~~~~~~~~~~~~

The AGORA model uses the :ref:`GEAR model <gear_star_formation>` scheme, however with the 
``GEARStarFormation:star_formation_mode`` parameter set to ``agora``. Instead of requiring the gas
density to reach the pressure floor, we simply require it to be denser than a density
threshold defined by ``GEARStarFormation:density_threshold``.


Recommended parameters for the AGORA model should be:


.. code:: YAML

  GEARStarFormation:
    star_formation_mode: agora            
    star_formation_efficiency: 0.01   
    maximal_temperature:  1e10       
    n_stars_per_particle: 1
    min_mass_frac: 0.5
    density_threshold:   1.67e-23   



.. _agora_feedback_and_chemistry:

Stellar Feedback and Chemistry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the current implementation, only two `elements`, iron (Fe) and the sum of all elements heavier than helium (metallicity)
are considered.
Only stellar feeback from core collapse supernovae (CCSNe) is considered.
The model assumes that a number ``AGORAFeedback:ccsne_per_solar_mass`` of CCSNe per solar mass will form out of 
each stellar particle. Those supernovae will all explode after a time 
``AGORAFeedback:supernovae_explosion_time_myr`` after the birth of the stellar particle.
At this time, each stellar particle will expel:

  - ``AGORAFeedback:ejected_mass_in_solar_mass_per_CCSN`` amount of gas (in solar mass), per CCSN formed
  - ``AGORAFeedback:ejected_Fe_mass_in_solar_mass_per_CCSN`` amount of iron (in solar mass), per CCSN formed
  - ``AGORAFeedback:ejected_metal_mass_in_solar_mass_per_CCSN`` amount of metals (in solar mass), per CCSN formed

In addition stellar particles will release energy:

  - ``AGORAFeedback:energy_in_erg_per_CCSN`` erg per CCSN formed

The energy released effectively into surrounding gas particles can be mitigated with the 
parameter ``AGORAFeedback:supernovae_efficiency``, used as a simple factor to ``AGORAFeedback:energy_in_erg_per_CCSN``.

Both energy and mass are ejected into the surrounding gas according to the SPH kernel.

The mass fraction of elements received by gas particles will be smoothed (smoothed metallicity scheme)
by using the SPH kernel. Note that this scheme does not exchange any material between particles. 
Snapshots can store both the smoothed metallicity (``SmoothedMetalMassFractions``)
and/or the non-smoothed one (``MetalMassFractions``), i.e., the mass fraction of elements effectively received
by the gas particles.


The initial metallicity of the gas can be defined by the parameter ``AGORAChemistry:initial_metallicity``.
A value less than 0 forces the code to take the gas metallicity from the initial condition file (snapshot).
Instead, if ``AGORAChemistry:scale_initial_metallicity`` is different than 0, the initial mass fraction
of elements will be set to:

  - ``AGORAChemistry:initial_metallicity`` time ``AGORAChemistry:solar_abundance_Metals`` for the iron
  - ``AGORAChemistry:initial_metallicity`` time ``AGORAChemistry:solar_abundance_Metals`` for the metals


Instead, they will be set to:

  - ``AGORAChemistry:initial_metallicity`` for the iron
  - ``AGORAChemistry:initial_metallicity`` for the metals



Recommended parameters for the AGORA model should be:

.. code:: YAML

   AGORAChemistry:
    initial_metallicity: 1  
    scale_initial_metallicity: 1            
    solar_abundance_Fe: 0.001771
    solar_abundance_Metals: 0.02  
    

.. code:: YAML
   
   AGORAFeedback:
    energy_in_erg_per_CCSN: 1e51
    supernovae_efficiency: 1
    supernovae_explosion_time_myr: 5
    ccsne_per_solar_mass : 0.010989
    ejected_mass_in_solar_mass_per_CCSN : 14.8
    ejected_Fe_mass_in_solar_mass_per_CCSN : 2.63
    ejected_metal_mass_in_solar_mass_per_CCSN : 2.63






.. _agora_pressure_floor:

Pressure Floor
~~~~~~~~~~~~~~

The AGORA model uses precisely the same pressure floor than the :ref:`GEAR model <gear_pressure_floor>`.



