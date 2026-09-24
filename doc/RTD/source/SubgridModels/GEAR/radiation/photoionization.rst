.. GEAR sub-grid model HII photoionization
   Darwin Roduit, 24th September 2026

.. _gear_radiation_hii:

HII photoionization
=====================

A young, hot star ionizes and heats the gas around it, carving out an HII region. GEAR models this as a Strömgren-sphere-like budget: each star is given an ionizing photon rate (:math:`Q_\mathrm{H}`, read from the radiation table), and it spends this budget on the gas particles around it, closest first.

A gas particle is eligible to be ionized if it is cold (:math:`T \leq 1.01 \times 10^4` K) and dense enough (physical density above ``GEARFeedback:HII_min_density_Hpcm3``). Claiming a fresh particle costs its remaining neutral hydrogen content, plus a maintenance reserve sized to the interval until the star's next pass; this reserve is what makes the region's growth implicit and numerically stable even when that interval is much longer than the local recombination time. Once claimed, a particle is held ionized as long as the star can afford its recombination losses every pass, computed from the case-B recombination coefficient (Hui and Gnedin 1997, a fit accurate to 0.7% from 1 K to :math:`10^9` K) evaluated at the temperature the ionized gas is held at. That temperature is the minimum of the energy needed to fully ionize the gas and a metallicity-dependent collisional-equilibrium temperature floor (Hopkins 2023's photoionization temperature fit). The floor itself falls from about :math:`6.6 \times 10^4` K at zero metallicity, down through :math:`10^4` K around :math:`Z \approx 0.23\,Z_\odot`, to lower values above solar metallicity, so a metal-poor particle is generally held hotter than a metal-rich one.

When a star's budget cannot fully afford its outermost (boundary) particle, ``GEARFeedback:HII_deterministic_boundary_ionization`` decides what happens: probabilistically (an unbiased coin flip, the default) or deterministically (always ionize it, letting the budget go slightly negative).

The star's search for gas to ionize is bounded by ``Stars:HII_max_search_radius`` (comoving), with automatic radius expansion and retry if a single pass does not reach every eligible particle. This search geometry is parsed whenever SWIFT is built with ``--with-feedback=GEAR``, and ``HII_max_search_radius`` has no default: a GEAR feedback run that never intends to use photoionization still has to set it.

Optionally, a star's ionizing budget can be split across several independent angular sectors (a HEALPix tessellation, ``GEARFeedback:HII_angular_nside``), so that one heavily-illuminated direction cannot starve another.

By default, an ionized particle is simply flagged and held at its floor temperature by the cooling module until its tag expires. Setting ``GEARFeedback:HII_couple_ionization_rate`` instead derives a physically motivated, distance-dependent photoionization rate coefficient (:math:`\Gamma_\mathrm{HI} = \sigma_\mathrm{HI} \times` ionizing flux, with :math:`\sigma_\mathrm{HI} = 6.3 \times 10^{-18}\,\mathrm{cm}^2` at the Lyman limit, Osterbrock and Ferland 2006) and feeds it into Grackle's own radiative-transfer rate fields every step, instead of the instant flag-and-floor scheme. This requires Grackle mode 1 or higher.

``GEARFeedback:HII_rebuild_time_Myr`` also bounds the timestep of every star still young enough to do photoionization (below ``HII_max_age_Myr``), not just the ionization cadence: a smaller value forces smaller, more numerous steps.

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

* ``HII_max_search_radius``: maximum comoving search radius a star extends to look for gas to ionize. Mandatory, no default.
* ``HII_max_retry_full_buffer``: number of full-buffer retries, at a fixed search radius, allowed per search pass. Default: 10.
* ``HII_max_radius_expansion_tries``: number of times the search radius is allowed to expand within one pass. Default: 5.
* ``HII_radius_expansion_factor``: growth factor applied to the search radius at each expansion try. Default: 1.1.

.. code:: YAML

   GEARFeedback:
     with_photoionization: 0                       # Master switch of the HII photoionization channel
     HII_min_density_Hpcm3: 1                       # Minimal density to consider a particle eligible for ionization (H atoms/cm^3)
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

Snapshot outputs
------------------

Two star fields are always written for a GEAR run with feedback, whether or not ``--with-tracers=GEAR`` is used:

.. list-table::
   :header-rows: 1
   :widths: 25 55 20

   * - Name
     - Description
     - Units
   * - ``HIIRegionRadii``
     - Comoving radius the star's HII region reached at its last budget rebuild
     - [U_L]
   * - ``HIIRegionMasses``
     - Gas mass the star currently holds ionized
     - [U_M]

These are the search algorithm's own bookkeeping, not a direct measurement of the gas's thermal state: previously-tagged gas can stay warm well past its tag's expiry without being re-tagged.

The following fields need ``--with-tracers=GEAR``:

.. list-table::
   :header-rows: 1
   :widths: 30 50 20

   * - Name
     - Description
     - Units
   * - ``IsIonizedFlags`` (gas)
     - Is this gas particle currently flagged as ionized?
     - [-]
   * - ``HIIStarIDs`` (gas)
     - ID of the star that ionized this particle
     - [-]
   * - ``FinalHIIRegionRadii`` (star)
     - ``HIIRegionRadii`` retired at the star's death or ineligibility
     - [U_L]
   * - ``FinalHIIRegionMasses`` (star)
     - ``HIIRegionMasses`` retired the same way
     - [U_M]

References
-----------

- Hui and Gnedin (1997), MNRAS 292, 27, Appendix A: case-B recombination coefficient.
- Hopkins (2023): metallicity-dependent photoionization temperature floor.
- Osterbrock and Ferland (2006): HI photoionization cross-section at the Lyman limit.
