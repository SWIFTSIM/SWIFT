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

.. warning::
   MPI support for this radiation model is work in progress. Test carefully before relying on a multi-rank run.

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
     - ``GEARFeedback:with_radiation_pressure`` (with ``GEARFeedback:radiation_pressure_efficiency``)
     - :ref:`gear_radiation_pressure`
   * - Interstellar radiation field
     - ``GEARFeedback:with_interstellar_radiation_field``
     - :ref:`gear_isrf`

Limitations
-------------

- **No sub-cycling.** Radiation rides the star's or gas particle's own timestep; there is no separate, finer radiation clock. See "Timestep criteria" below.
- **Two ISRF bands only.** The ISRF module transports the PE (6.0 to 11.2 eV) and Lyman-Werner (11.2 to 13.6 eV) bands. Radiation pressure and HII photoionization instead use the star's bolometric luminosity and ionizing photon rate as single lumped quantities, not a spectrum: no channel resolves the field by wavelength outside the two ISRF bands.
- **H2 photodissociation needs Grackle mode 2 or 3.** Below that, the Lyman-Werner band is still computed, transported and written to snapshots, and it still contributes to photoelectric heating, but it does not dissociate :math:`\mathrm{H}_2`. See the Grackle cooling mode paragraph below.
- **HII rate-coupling is H-only.** ``GEARFeedback:HII_couple_ionization_rate`` derives a per-particle HI photoionization rate; the HeI/HeII ionization rates Grackle also accepts stay at their separate, spatially-uniform ``GrackleCooling`` scalars.
- **No ray tracing or shadowing.** HII photoionization spends a star's photon budget on gas by distance (or, with ``HII_angular_nside``, by angular sector), not along a traced sightline: intervening dense gas does not shield a farther particle. Without ``ISRF_propagation``, the injected ISRF field only reaches gas inside the illuminating star's own SPH kernel; with it on, the field is transported by a hyperbolic moment (M1) scheme at a reduced speed of light, not a ray trace.

Configuring and compiling
--------------------------

The radiation model is part of the GEAR feedback module, so it is built whenever GEAR feedback is:

.. code:: bash

   ./configure --with-chemistry=GEAR_10 --with-feedback=GEAR \
               --with-cooling=grackle_0 --with-stars=GEAR --with-sink=GEAR \
               --with-star-formation=GEAR --with-tracers=GEAR \
               --with-kernel=wendland-C2 --with-grackle=path/to/grackle

Radiation itself only needs ``--with-feedback=GEAR``, its matching ``--with-stars=GEAR`` (the mandatory ``Stars:HII_max_search_radius`` parameter and the HII search-radius machinery exist only in the GEAR stars particle, so no other stars module will compile against it), and a Grackle cooling mode (``--with-cooling=grackle_N``, ``N`` chosen per the requirements below and in :ref:`gear_isrf`). ``--with-chemistry=GEAR_10`` is not itself required to build the model, but every metallicity-dependent term on this page (the HII temperature floor, the radiation-pressure opacities, the ISRF dust extinction) reads the particle's metallicity, so a non-GEAR chemistry model gives every particle ``Z=0`` and silently disables those terms. ``--with-sink=GEAR``, ``--with-star-formation=GEAR`` and ``--with-kernel=wendland-C2`` are not read by the radiation code at all: they come from the shipped examples' own full-simulation configure line, needed to form and evolve the stars that radiation then acts on, not required by the radiation channels themselves. ``--with-tracers=GEAR`` is read in exactly one place: radiation pressure's own momentum and kick diagnostics (described next). No other radiation channel touches the tracers module. This is why :ref:`gear_isrf`'s own configure line lists only five options, marked "at least": it is the radiation-only minimum, not a full production build.

**The GEAR tracers module is optional to build.** Every radiation channel compiles without ``--with-tracers=GEAR``: the HII ionization tag lives outside the tracers module, and none of the HII or ISRF snapshot fields are tracers-module fields, so none of them need ``--with-tracers`` (see :ref:`gear_output_hii` and :ref:`gear_output_isrf`). The one exception is radiation pressure's own diagnostic fields, ``CumulativeMomentumFromRadiationPressure`` and ``MaxKickVelocityFromRadiationPressure`` (see :ref:`gear_output_radiation_pressure`): those still live in the GEAR tracers module, so build with ``--with-tracers=GEAR`` if you want them. Every GEAR feedback build still requires the mandatory ``Stars:HII_max_search_radius`` parameter in the YAML file (see :ref:`gear_radiation_hii`), even for a run with photoionization switched off.

**The Grackle cooling mode gates the ISRF's two effects differently.** ``--with-cooling=grackle_N`` sets how many chemical species Grackle tracks (0: none, 1: H/He, 2: + :math:`\mathrm{H}_2`, 3: + D). HII photoionization and radiation pressure work at any mode. For the ISRF:

- photoelectric heating and dust chemistry are fed by the sum of the PE and Lyman-Werner band energies (the Habing field), at any Grackle mode, ``grackle_0`` included;
- :math:`\mathrm{H}_2` photodissociation is fed by the Lyman-Werner band alone, and only at ``grackle_2`` or ``grackle_3``, since :math:`\mathrm{H}_2` is untracked below that.

So at ``grackle_0`` or ``grackle_1``, the Lyman-Werner band is still computed, transported and written to snapshots, and its energy still counts towards photoelectric heating through the Habing sum above. What is inert is only its own dedicated effect: it never dissociates :math:`\mathrm{H}_2` below mode 2, and SWIFT gives no start-up warning for this case. Rate-coupling the HII ionizing rate into Grackle (``GEARFeedback:HII_couple_ionization_rate``) requires ``grackle_1`` or higher.

**The** ``--with-subgrid=GEAR`` **and** ``--with-subgrid=GEAR-G3`` **shortcuts fix both the Grackle mode and the Jeans pressure floor** (:ref:`gear_pressure_floor`) as a pair, and configure rejects either shortcut combined with an explicit ``--with-cooling`` or ``--with-pressure-floor`` option. ``--with-subgrid=GEAR`` gives ``grackle_0`` with the pressure floor on: it can never dissociate :math:`\mathrm{H}_2` through the Lyman-Werner band, with no warning. ``--with-subgrid=GEAR-G3`` gives ``grackle_3`` with the pressure floor off. For both :math:`\mathrm{H}_2` photodissociation and the pressure floor, drop ``--with-subgrid`` and pass the individual options shown above, with ``--with-cooling=grackle_2`` or ``grackle_3`` plus ``--with-pressure-floor=GEAR`` (the individual ``--with-pressure-floor`` option defaults to ``none`` if left out, same as ``GEAR-G3``).

Running an example
-------------------

Each channel has its own worked examples under ``examples/SubgridTests/StellarFeedback/``, each with a README covering its own configure line, run command and check script. There is no radiation-specific command-line flag: SWIFT decides which channels run from the ``GEARFeedback:`` parameters in the table above, once the run itself has stars and feedback switched on. A minimal invocation, once the parameter file sets the channels you want, is:

.. code:: bash

   ./swift --hydro --self-gravity --stars --star-formation --feedback --cooling --threads=4 params.yml

The ``--gear`` shortcut expands to exactly ``--hydro --limiter --sync --self-gravity --stars --star-formation --cooling --feedback``; note that it does not add ``--sinks``, so a run that also wants sink particles needs that flag on top. Each shipped example's own ``run.sh`` adjusts the flag set to its own setup (some also drop gravity in favour of ``--external-gravity``, or leave out cooling to isolate one channel); see the example's README for the exact flags it uses.

- ``HIIRegions/Starbench`` is a good place to start with photoionization: a single ionizing source, checked against the STARBENCH D-type expansion solution of Bisbas et al. (2015). ``HIIRegions/StromgrenSphere`` is a simpler setup with the same geometry, but its own README states that it exists to exercise the search-radius/rebuild-cadence machinery rather than to reproduce a published result.
- ``RadiationPressure/RadiationPressureShellExpansion`` exercises the radiation-pressure channel on an expanding shell, checked against the analytic solution of Krumholz and Matzner (2009).
- ``ISRF/ISRFPhotoelectricHeating`` is the simplest ISRF setup: photoelectric heating from a single star, no propagation. If you specifically want to check :math:`\mathrm{H}_2` photodissociation, use ``ISRF/ISRFH2Photodissociation`` instead.

The full set under ``HIIRegions/`` and ``ISRF/`` covers more specific behaviour (cosmology, propagation schemes, dissipation), one aspect per directory.

Timestep criteria
--------------------

A star still young enough to do photoionization (below ``HII_max_age_Myr``) has its timestep bounded by ``GEARFeedback:HII_rebuild_time_Myr``, so it rebuilds its HII region on schedule; see :ref:`gear_radiation_hii`. Every star, whatever channel it runs, also has its event-anchored timestep terms floored by ``GEARFeedback:event_dt_floor_Myr``, and a Single Stellar Population or continuous-IMF star additionally has its own timestep tightened while young by ``GEARFeedback:dt_evolution_factor_max``; see :ref:`gear_stellar_evolution_and_feedback`.

On the gas side, a particle generally rides its own hydrodynamic timestep: no radiation channel adds a general CFL-like bound to it. The one exception is the ISRF's fixed-fraction-of-c propagation scheme (``ISRF_c_hyp_scheme: 2``), where ``GEARFeedback:ISRF_c_hyp_fixed_fraction_of_c`` adds a receiver-side stability bound to a particle carrying or near the field; see :ref:`gear_isrf`. A gas particle freshly ionized or freshly illuminated is synced onto a shorter time bin on the step it is first touched, but this is a one-off wake-up, not sub-cycling: it then rides its own step like any other particle.

Checking a run
----------------

At start-up SWIFT logs each channel's on/off state, and, when a channel is on, its main efficiency or margin value (grep the log for ``Photoionization``, ``Radiation pressure`` or ``Photo-electric heating``), so a switch left off by mistake shows up before the run gets far. To stop a production snapshot from carrying every radiation diagnostic field, select the fields you actually want with an output-selection file; see :ref:`Output_selection_label`. The full list of radiation snapshot fields is on the :ref:`gear_output_fields` page.
