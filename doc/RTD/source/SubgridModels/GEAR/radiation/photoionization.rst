.. GEAR sub-grid model HII photoionization
   Darwin Roduit, 24th September 2026

.. _gear_radiation_hii:

HII photoionization
=====================

A young, hot star ionizes and heats the gas around it, carving out an HII region. GEAR models this as a Strömgren (Stromgren)-sphere-like budget: each star is given an ionizing photon rate (:math:`Q_\mathrm{H}`, read from the radiation table), and it spends this budget on the gas particles around it, closest first.

A gas particle is eligible to be ionized if it is cold (:math:`T \leq 1.01 \times 10^4` K) and dense enough (physical density above ``GEARFeedback:HII_min_density_Hpcm3``). Once ionized, a particle is held at a temperature that is the minimum of the energy needed to fully ionize the gas and a metallicity-dependent collisional-equilibrium temperature floor (Hopkins 2023's photoionization temperature fit). The floor itself falls from about :math:`6.6 \times 10^4` K at zero metallicity, down through :math:`10^4` K around :math:`Z \approx 0.23\,Z_\odot`, to lower values above solar metallicity, so a metal-poor particle is generally held hotter than a metal-rich one.

Claiming a fresh particle costs its remaining neutral hydrogen content, plus a maintenance reserve sized to the interval until the star's next pass; this reserve is what makes the region's growth implicit and numerically stable even when that interval is much longer than the local recombination time. A particle stays ionized as long as the star can afford its recombination losses every pass, computed from the case-B recombination coefficient (Hui and Gnedin 1997, a fit accurate to 0.7% from 1 K to :math:`10^9` K) evaluated at the held temperature above.

When a star's budget cannot fully afford its outermost (boundary) particle, ``GEARFeedback:HII_deterministic_boundary_ionization`` decides what happens: probabilistically (an unbiased coin flip, the default) or deterministically (always ionize it, letting the budget go slightly negative).

The star's search for gas to ionize is bounded by ``Stars:HII_max_search_radius`` (comoving), with automatic radius expansion and retry if a single pass does not reach every eligible particle. This search geometry is parsed whenever SWIFT is built with ``--with-feedback=GEAR``, and ``HII_max_search_radius`` has no default: a GEAR feedback run that never intends to use photoionization still has to set it.

Optionally, a star's ionizing budget can be split across several independent angular sectors (a HEALPix tessellation, ``GEARFeedback:HII_angular_nside``), so that one heavily-illuminated direction cannot starve another.

By default, an ionized particle is simply flagged and held at its floor temperature by the cooling module until its tag expires. Setting ``GEARFeedback:HII_couple_ionization_rate`` instead derives a physically motivated, distance-dependent photoionization rate coefficient (:math:`\Gamma_\mathrm{HI} = \sigma_\mathrm{HI} \times` ionizing flux, with :math:`\sigma_\mathrm{HI} = 6.3 \times 10^{-18}\,\mathrm{cm}^2` at the Lyman limit, Osterbrock and Ferland 2006) and feeds it into Grackle's own radiative-transfer rate fields every step, instead of the instant flag-and-floor scheme. This requires Grackle mode 1 or higher.

``GEARFeedback:HII_rebuild_time_Myr`` also bounds the timestep of every star still young enough to do photoionization (below ``HII_max_age_Myr``), not just the ionization cadence: a smaller value forces smaller, more numerous steps.

Configuring and compiling
--------------------------

See :ref:`gear_radiation` for how to configure and build this channel, including the ``--with-tracers=GEAR`` requirement shared by every radiation channel.

Model parameters
------------------

* ``with_photoionization``: master switch. 1 turns the channel on, 0 (default) leaves it off.
* ``HII_min_density_Hpcm3``: minimum physical gas density (hydrogen atoms/cm^3) for a particle to be eligible for ionization. Default: 1. The default suits dense star-forming gas; a lower-density idealised box (a uniform test box, for instance) needs a lower value, or nothing is ever ionized. The shipped ``StromgrenSphere``, ``TaskGraphPair`` and ``TaskGraphCluster`` examples set it to 0.01.
* ``HII_max_age_Myr``: age above which a star stops doing photoionization. Default: 50.
* ``HII_rebuild_time_Myr``: interval between HII budget rebuilds. A negative value rebuilds every step instead of on a fixed schedule. Default: 0.5.
* ``HII_rebuild_floor_Myr``: floor on the interval the per-pass photon budget is integrated over. Must be positive. Default: 1e-4.
* ``HII_angular_nside``: HEALPix splitting of the ionizing budget. 0 (default) is spherical, one shared budget; N >= 1 gives :math:`12 N^2` angular pixels. N > 0 requires building with ``--with-chealpix``, and the resulting pixel count must fit the build's ``--with-number-of-hii-angular-pixels`` (default 12, i.e. N <= 1).
* ``HII_deterministic_boundary_ionization``: 0 (default, probabilistic) or 1 (deterministic) treatment of the boundary particle described above.
* ``HII_couple_ionization_rate``: rate-couple the photoionization rate into Grackle instead of the default flag-and-floor scheme. 0 (default) or 1. Requires ``--with-cooling=grackle_1`` or higher; forces ``GrackleCooling:use_radiative_transfer`` on internally.

Four further parameters, in the ``Stars`` section, control the search geometry and are parsed whenever SWIFT is built with GEAR feedback:

* ``HII_max_search_radius``: maximum search radius a star extends to look for gas to ionize, in comoving internal units (physical reach = ``a`` times this value). Mandatory, no default. This is a hard cap, not just a starting value for the expansion below: if the true Strömgren radius would need a larger search, the region simply stops growing at this radius, and gas beyond it is never ionized, with no error. Raise it if you expect a physically larger HII region than the cap allows.
* ``HII_max_retry_full_buffer``: number of full-buffer retries, at a fixed search radius, allowed per search pass. Default: 10.
* ``HII_max_radius_expansion_tries``: number of times the search radius is allowed to expand within one pass. Default: 5.
* ``HII_radius_expansion_factor``: growth factor applied to the search radius at each expansion try. Default: 1.1.

.. code:: YAML

   GEARFeedback:
     with_photoionization: 0                       # Master switch of the HII photoionization channel
     HII_min_density_Hpcm3: 1                       # Minimal density to consider a particle eligible for ionization (H atoms/cm^3). Lower this for a low-density idealised box, or nothing is ever ionized
     HII_max_age_Myr: 50                            # Age (Myr) at which stars stop doing photoionization
     HII_rebuild_time_Myr: 0.5                      # Time (Myr) between HII budget rebuilds; negative rebuilds every step
     HII_rebuild_floor_Myr: 1e-4                    # Floor (Myr) on the interval the budget is integrated over
     HII_angular_nside: 0                           # HEALPix splitting of the ionizing budget: 0 = spherical, N>=1 = 12*N^2 pixels
     HII_deterministic_boundary_ionization: 0       # 0 = probabilistic, 1 = deterministic boundary-particle treatment
     HII_couple_ionization_rate: 0                  # Rate-couple into Grackle's own RT fields instead of flag-and-floor

   Stars:
     HII_max_search_radius: 0.1                     # Maximal comoving search radius for gas to ionize
     HII_max_retry_full_buffer: 10                  # (Optional) Full-buffer retries per search pass (Default: 10)
     HII_max_radius_expansion_tries: 5               # (Optional) Search-radius expansion tries per pass (Default: 5)
     HII_radius_expansion_factor: 1.1                # (Optional) Growth factor per expansion try (Default: 1.1)

Compile-time constants
------------------------

``--with-number-of-hii-angular-pixels`` (default 12) caps the HEALPix pixel count ``HII_angular_nside`` may request; see above.

Three ``-D`` flags, set via ``CFLAGS+="..." ./configure``, are for test and diagnostic setups, not production runs. The shipped ``Starbench``, ``ClumpSmith2021``, ``Hu2017`` and ``StromgrenSphereCosmo`` examples use them to build the idealised, two-state configuration their check scripts compare against:

* ``-DIONIZATION_FEEDBACK_DEBUG_NO_COOLING``: tagged gas never cools, so a tag never expires.
* ``-DIONIZATION_FEEDBACK_DEBUG_FIXED_IONIZED_TEMPERATURE_K=<K>``: holds every ionized particle at a fixed temperature regardless of metallicity, instead of the metallicity-dependent floor above.
* ``-DIONIZATION_FEEDBACK_DEBUG_FIXED_NEUTRAL_TEMPERATURE_K=<K>``: also pins neutral gas to a fixed temperature. Combined with the previous flag, this gives an idealised two-temperature medium.

Snapshot outputs
------------------

The gas and star fields this channel writes are documented on the :ref:`gear_output_hii` section of the :ref:`gear_output_fields` page.

References
-----------

- Hui and Gnedin (1997), MNRAS 292, 27, Appendix A: case-B recombination coefficient.
- Hopkins (2023): metallicity-dependent photoionization temperature floor.
- Osterbrock and Ferland (2006): HI photoionization cross-section at the Lyman limit.
