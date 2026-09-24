.. GEAR sub-grid model radiation feedback
   Darwin Roduit, 24th September 2026

.. _gear_radiation:

Radiation feedback
===================

Young, massive stars couple to the surrounding gas through radiation long before their first supernova. GEAR models three separate channels of this coupling:

- **HII photoionization**: a star ionizes and heats the gas around it, mimicking a Strömgren sphere.
- **Radiation pressure**: the star's bolometric light pushes on the dusty gas around it through momentum deposition.
- **Interstellar radiation field (ISRF)**: a non-ionizing ultraviolet field that heats dust (photoelectric heating) and dissociates :math:`\mathrm{H}_2` (Lyman-Werner), optionally transported away from the source with a hyperbolic moment scheme.

Every channel reads the star's photon output (an ionizing photon rate, a bolometric luminosity, band luminosities) from the ``Data/Radiation`` group of the stellar evolution table (``GEARFeedback:yields_table``), as a function of the star's mass, metallicity and age. See :ref:`gear_radiation_tables` for what this group must contain and how to get one. The three channels are otherwise independent: each is switched on by its own parameter, and you can run any combination of them.

- :ref:`gear_radiation_hii`
- :ref:`gear_radiation_pressure`
- :ref:`gear_isrf`
- :ref:`gear_radiation_tables`

Configuring and compiling
--------------------------

The radiation model is part of the GEAR feedback module, so it is built whenever GEAR feedback is:

.. code:: bash

   ./configure --with-chemistry=GEAR_10 --with-feedback=GEAR \
               --with-cooling=grackle_0 --with-stars=GEAR --with-sink=GEAR \
               --with-star-formation=GEAR --with-tracers=GEAR \
               --with-kernel=wendland-C2 --with-grackle=path/to/grackle

**--with-tracers=GEAR is required, not optional.** The HII ionization tag a gas particle carries lives in the GEAR tracers module's own per-particle data. SWIFT will not compile ``--with-feedback=GEAR`` without ``--with-tracers=GEAR`` alongside it, whatever radiation channel you actually intend to use.

**The Grackle cooling mode limits what you get.** ``--with-cooling=grackle_N`` sets how many chemical species Grackle tracks (0: none, 1: H/He, 2: + :math:`\mathrm{H}_2`, 3: + D). HII photoionization and radiation pressure work at any mode. The ISRF's Lyman-Werner channel needs :math:`\mathrm{H}_2` to dissociate, so it requires ``grackle_2`` or ``grackle_3``; at a lower mode the photoelectric-heating half of the ISRF still runs. Rate-coupling the HII ionizing rate into Grackle (``GEARFeedback:HII_couple_ionization_rate``) requires ``grackle_1`` or higher.

Running an example
-------------------

Each channel has its own worked examples under ``examples/SubgridTests/StellarFeedback/``, each with a README covering its own configure line, run command and check script:

- ``HIIRegions/StromgrenSphere`` is the simplest place to start with photoionization: a single star in a uniform box, checked against the classical Spitzer solution.
- ``RadiationPressure/RadiationPressureShellExpansion`` exercises the radiation-pressure channel on an expanding shell.
- ``ISRF/ISRFPhotoelectricHeating`` is the simplest ISRF setup: photoelectric heating from a single star, no propagation.

The full set under ``HIIRegions/`` and ``ISRF/`` covers more specific behaviour (cosmology, task-graph cadence, propagation schemes, dissipation), one aspect per directory.
