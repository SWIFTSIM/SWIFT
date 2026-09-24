.. GEAR sub-grid model radiation feedback
   Darwin Roduit, 24th September 2026

.. _gear_radiation:

Radiation feedback
===================

Young, massive stars couple to the surrounding gas through radiation long before their first supernova. GEAR models three separate channels of this coupling:

- **HII photoionization**: a star ionizes and heats the gas around it, mimicking a Strömgren (Stromgren) sphere.
- **Radiation pressure**: the star's bolometric light pushes on the dusty gas around it through momentum deposition.
- **Interstellar radiation field (ISRF)**: a non-ionizing ultraviolet field that heats dust (photoelectric heating) and dissociates :math:`\mathrm{H}_2` (Lyman-Werner), optionally transported away from the source with a hyperbolic moment scheme.

This is a sub-grid coupling built into the GEAR feedback module. It is a different code path from the mesh-based radiative-transfer solver at :ref:`rt_GEAR` (``--with-rt=GEAR_N``): that solver replaces the hydro scheme with GIZMO-MFV and solves photon transport on the mesh, while the model on this page runs inside any ordinary GEAR feedback build and only moves energy and ionization state between a star and its SPH neighbours.

Every channel reads the star's photon output (an ionizing photon rate, a bolometric luminosity, band luminosities) from the ``Data/Radiation`` group of the stellar evolution table (``GEARFeedback:yields_table``), as a function of the star's mass, metallicity and age. See :ref:`gear_radiation_tables` for what this group must contain and how to get one. The three channels are otherwise independent: each is switched on by its own parameter, and you can run any combination of them.

.. list-table:: Switching a channel on
   :header-rows: 1
   :widths: 30 40 30

   * - Channel
     - Switch
     - Page
   * - HII photoionization
     - ``GEARFeedback:with_photoionization``
     - :ref:`gear_radiation_hii`
   * - Radiation pressure
     - ``GEARFeedback:radiation_pressure_efficiency`` (0 = off)
     - :ref:`gear_radiation_pressure`
   * - Interstellar radiation field
     - ``GEARFeedback:with_interstellar_radiation_field``
     - :ref:`gear_isrf`

Configuring and compiling
--------------------------

The radiation model is part of the GEAR feedback module, so it is built whenever GEAR feedback is:

.. code:: bash

   ./configure --with-chemistry=GEAR_10 --with-feedback=GEAR \
               --with-cooling=grackle_0 --with-stars=GEAR --with-sink=GEAR \
               --with-star-formation=GEAR --with-tracers=GEAR \
               --with-kernel=wendland-C2 --with-grackle=path/to/grackle

**The convenience option** ``--with-subgrid=GEAR`` **builds this too, but it silently pins the cooling mode to** ``grackle_0``, **which drops the Lyman-Werner band** (see the Grackle cooling mode paragraph below). If you want the Lyman-Werner channel without configuring every option by hand, use ``--with-subgrid=GEAR-G3`` instead: it selects ``grackle_3``. Configuring the options individually, as in the command above, lets you pick any Grackle mode directly.

Radiation itself only needs ``--with-feedback=GEAR``, its matching ``--with-stars=GEAR`` (the mandatory ``Stars:HII_max_search_radius`` parameter and the HII search-radius machinery exist only in the GEAR stars particle, so no other stars module will compile against it), ``--with-tracers=GEAR``, and a Grackle cooling mode (``--with-cooling=grackle_N``, ``N`` chosen per the requirements below and in :ref:`gear_isrf`). ``--with-chemistry=GEAR_10`` is not itself required to build the model, but every metallicity-dependent term on this page (the HII temperature floor, the radiation-pressure opacities, the ISRF dust extinction) reads the particle's metallicity, so a non-GEAR chemistry model gives every particle ``Z=0`` and silently disables those terms. ``--with-sink=GEAR``, ``--with-star-formation=GEAR`` and ``--with-kernel=wendland-C2`` are not read by the radiation code at all: they come from the shipped examples' own full-simulation configure line, needed to form and evolve the stars that radiation then acts on, not required by the radiation channels themselves. This is why :ref:`gear_isrf`'s own configure line lists only six options, marked "at least": it is the radiation-only minimum, not a full production build.

**The GEAR tracers module is required, not optional.** The HII ionization tag a gas particle carries lives in the GEAR tracers module's own per-particle data. SWIFT will not compile ``--with-feedback=GEAR`` without ``--with-tracers=GEAR`` alongside it, whatever radiation channel you actually intend to use, and every GEAR feedback build also requires the mandatory ``Stars:HII_max_search_radius`` parameter in the YAML file (see :ref:`gear_radiation_hii`), even for a run with photoionization switched off.

**The Grackle cooling mode limits what you get.** ``--with-cooling=grackle_N`` sets how many chemical species Grackle tracks (0: none, 1: H/He, 2: + :math:`\mathrm{H}_2`, 3: + D). HII photoionization and radiation pressure work at any mode. The ISRF's Lyman-Werner channel needs :math:`\mathrm{H}_2` to dissociate, so it requires ``grackle_2`` or ``grackle_3``; at a lower mode the photoelectric-heating half of the ISRF still runs. Rate-coupling the HII ionizing rate into Grackle (``GEARFeedback:HII_couple_ionization_rate``) requires ``grackle_1`` or higher.

Running an example
-------------------

Each channel has its own worked examples under ``examples/SubgridTests/StellarFeedback/``, each with a README covering its own configure line, run command and check script. There is no radiation-specific command-line flag: SWIFT decides which channels run from the ``GEARFeedback:`` parameters in the table above, once the run itself has stars and feedback switched on. A minimal invocation, once the parameter file sets the channels you want, is:

.. code:: bash

   ./swift --hydro --self-gravity --stars --star-formation --feedback --cooling --threads=4 params.yml

The ``--gear`` shortcut expands to exactly ``--hydro --limiter --sync --self-gravity --stars --star-formation --cooling --feedback``; note that it does not add ``--sinks``, so a run that also wants sink particles needs that flag on top. Each shipped example's own ``run.sh`` adjusts the flag set to its own setup (some also drop gravity in favour of ``--external-gravity``, or leave out cooling to isolate one channel); see the example's README for the exact flags it uses.

- ``HIIRegions/Starbench`` is a validated place to start with photoionization: a single ionizing source, checked against the STARBENCH D-type expansion solution of Bisbas et al. (2015). ``HIIRegions/StromgrenSphere`` is a simpler setup with the same geometry, but its own README states that it exists to exercise the search-radius/rebuild-cadence machinery rather than to reproduce a published result.
- ``RadiationPressure/RadiationPressureShellExpansion`` exercises the radiation-pressure channel on an expanding shell, checked against the analytic solution of Krumholz and Matzner (2009).
- ``ISRF/ISRFPhotoelectricHeating`` is the simplest ISRF setup: photoelectric heating from a single star, no propagation. If you specifically want to check :math:`\mathrm{H}_2` photodissociation, use ``ISRF/ISRFH2Photodissociation`` instead.

The full set under ``HIIRegions/`` and ``ISRF/`` covers more specific behaviour (cosmology, propagation schemes, dissipation), one aspect per directory.

Checking a run
----------------

At start-up SWIFT logs each channel's on/off state, and, when a channel is on, its main efficiency or margin value (grep the log for ``Photoionization``, ``Radiation pressure`` or ``Photo-electric heating``), so a switch left off by mistake shows up before the run gets far. To stop a production snapshot from carrying every radiation diagnostic field, select the fields you actually want with an output-selection file; see :ref:`Output_selection_label`.
