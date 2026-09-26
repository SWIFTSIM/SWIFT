.. GEAR stellar evolution and radiation table
   Darwin Roduit, 25th September 2026

.. _gear_stellar_evolution_table:

The stellar evolution table
============================

GEAR reads almost every stellar-evolution and radiation quantity from a single HDF5 file, not from the YAML parameter file. The file is reached by ``GEARFeedback:yields_table`` (population II stars) and ``GEARFeedback:yields_table_first_stars`` (population III stars, used below the metallicity set by ``GEARFeedback:imf_transition_metallicity``). This page describes what the file must contain and how to get or check one. To generate your own table, use pychem's ``pychem_generate_hdf5_parameters``; refer to the pychem documentation for the parameter file itself.

Table structure
----------------

.. graphviz:: feedback_table.dot

Solid (dashed) squares represent a group (a dataset), with the object's name underlined and its attributes written below. Everything in the groups below is in solar masses or is unitless (a mass fraction or an exponent); the ``Data/Radiation`` group, described in its own section below, is the exception and carries physical CGS units instead.

* ``Data``: ``elts`` is the array of element names (the last entry, ``Metals``, is the sum of every other element), ``MeanWDMass`` is the white dwarf mass, and ``SolarMassAbundances`` is the solar mass fraction of each element.
* ``IMF``: ``n + 1`` is the number of IMF segments, ``as`` their exponents (``n + 1`` values), ``ms`` the mass limits between segments (``n`` values), and ``Mmin``/``Mmax`` the minimal/maximal stellar mass. GEAR uses the Kroupa (2001) IMF restricted to ``[Mmin, Mmax]``.
* ``LifeTimes``: the stellar lifetime fit's coefficients, as a single 3x3 ``coeff_z`` table (Poirier 2004, summarised in Hausammann 2021).
* ``SNIa``: ``a`` is the binary-distribution exponent, ``bb1``/``bb2`` the companion probabilities :math:`b_i`, and the remaining attributes follow the SNIa rate formula's own names. Its ``Metals`` subgroup carries the element names (``elts``, matching ``Data``) and the metal mass fraction ejected per supernova (``data``).
* ``SNII``: mass limits ``Mmin``/``Mmax``. The yields are ``Ej`` (processed ejected mass fraction), ``Ejnp`` (non-processed) and one dataset per element in ``elts``, each uniformly sampled in log-mass with attributes ``min`` (first mass, in log) and ``step``.
* ``SW`` (stellar winds, Deng et al. 2024b): under ``MetallicityDependent``, four 2D datasets, ``Energy``/``Mass`` (per single star) and ``Integrated_Energy``/``Integrated_Mass_Loss`` (per SSP or continuous-star particle), each on a mass x metallicity grid described by its own ``dims``, ``m0``/``z0``, ``dm``/``dz`` and ``nm``/``nz`` attributes.

See :ref:`gear_stellar_evolution_and_feedback` for the physics these groups drive (the IMF, lifetimes, SNII/SNIa rates and yields, stellar winds).

.. _gear_radiation_tables:

The ``Data/Radiation`` group
------------------------------

Photoionization, radiation pressure and the interstellar radiation field (ISRF) read the star's photon output from a ``Data/Radiation`` group in the same file, as a function of the star's mass (a mass-only, "M" table) or mass and metallicity (a mass x metallicity, "M,Z" table). A table without this group cannot drive any radiation channel: SWIFT stops at start-up. Every dataset carries its own CGS ``units`` attribute, which SWIFT converts to internal units on read; a dataset without one is rejected.

Whenever any radiation channel is switched on, the group must carry:

* ``Luminosity`` and ``Integrated_Luminosity`` (erg/s and erg/s/Msun of stars formed): bolometric luminosity, used by radiation pressure.
* ``Q_H`` and ``Integrated_Q_H`` (1/s and 1/s/Msun): ionizing photon rate, used by HII photoionization.
* ``DotEExcess`` and ``Integrated_DotEExcess`` (erg/s and erg/s/Msun): excess-photon-energy emission rate above the HI ionization threshold.
* ``MeanPhotonEnergyLW`` and ``Integrated_MeanPhotonEnergyLW`` (both erg): photon-number-weighted mean Lyman-Werner photon energy. Required unconditionally, even for a photoionization-only or radiation-pressure-only run that never touches the ISRF. Unlike the other ``Integrated_`` datasets, this one is not a per-Msun rate: it is an intensive per-photon quantity (:math:`L_\mathrm{LW}/Q_\mathrm{LW}`), so it does not scale with the mass formed.

If the ISRF module is on (``GEARFeedback:with_interstellar_radiation_field: 1``), the group must also carry ``L_PE``/``Integrated_L_PE`` and ``L_LW``/``Integrated_L_LW`` (erg/s and erg/s/Msun), the photoelectric and Lyman-Werner band luminosities.

``Teff`` (photospheric effective temperature, K) is optional and has no ``Integrated_`` counterpart; it is read if present and otherwise simply unavailable.

Getting a table
------------------

``examples/GEAR_ICs_and_SCRIPTS/getChemistryTable.sh`` fetches the base table (stellar evolution only). By itself, or with ``--with-winds``, it downloads a table whose ``Data/Radiation`` group (if any) predates ``MeanPhotonEnergyLW``: SWIFT rejects it for every radiation channel, photoionization included. Add ``--with-radiation`` to also fetch the tables below that do carry usable radiation data, or run ``getRadiationTable.sh`` directly:

.. list-table:: Published tables
   :header-rows: 1
   :widths: 25 15 20 40

   * - Table
     - Axis
     - ``Data/Radiation``
     - Notes
   * - ``PopII_parsec_spectral.hdf5``
     - M,Z
     - Full (spectral, including ``L_PE``/``L_LW``)
     - pychem/PARSEC spectral synthesis. **Use this one for a science run.** All shipped ISRF examples use it.
   * - ``PopIII_parsec_spectral.hdf5``
     - M,Z
     - Full (spectral)
     - Same as above, for population III (first) stars.
   * - ``radiation_fits_popII.hdf5``
     - M
     - Full, but blackbody-derived, not spectral
     - Mass-only, blackbody-fits table. Its ``L_PE``/``L_LW`` come from integrating a blackbody spectrum over each band rather than a spectral synthesis library. The shipped HII-region and radiation-pressure examples require this exact file: their star masses were chosen to hit this table's own ``Q_H`` at a specific value (for instance, Starbench's 26.75 Msun reproduces Bisbas et al. (2015)'s :math:`10^{49}` photons/s), so swapping in the spectral table would change those examples' checked numbers.
   * - ``POPIIsw.h5`` / ``POPII.hdf5`` / ``POPIII_PISNe.hdf5``
     - --
     - Missing or incomplete
     - Fetched by ``getChemistryTable.sh`` without ``--with-radiation``. ``POPII.hdf5``/``POPIII_PISNe.hdf5`` (``--with-winds``) carry no ``Data/Radiation`` group at all. ``POPIIsw.h5`` (default) carries one, but it predates ``MeanPhotonEnergyLW``/``Integrated_MeanPhotonEnergyLW``, so it still fails the check above. Neither can drive any radiation channel.

Checking a table
------------------

``examples/GEAR_ICs_and_SCRIPTS/checkRadiationTable.sh <table.h5> [--with-isrf] [--require-1d]`` verifies a table before you spend time on the rest of an example's setup:

* with no flag, it checks the datasets every radiation channel needs (``Luminosity``, ``Q_H``, ``DotEExcess``, ``MeanPhotonEnergyLW`` and their ``Integrated_`` counterparts);
* ``--with-isrf`` additionally checks the four ISRF band datasets (``L_PE``, ``L_LW`` and their ``Integrated_`` counterparts);
* ``--require-1d`` additionally fails on a mass x metallicity ("M,Z") table, for an example whose own Python check script only understands a mass-only table.

It does not check the stellar-evolution groups above (``IMF``, ``LifeTimes``, ``SNII``, ``SNIa``, ``SW``); those are validated by SWIFT itself reading the table at start-up.

Interpolation
---------------

``interpolation_size_mass`` (section ``GEARRadiation``): number of points in the mass-axis interpolation grid built from the ``Data/Radiation`` table. Must be at least 2. Default: 500, matching the shipped tables' own native grid size; a smaller value costs interpolation error near the source fits' own mass-relation breakpoints.

.. code:: YAML

   GEARRadiation:
     interpolation_size_mass: 500   # Number of points in the mass interpolation of the Data/Radiation table (Default: 500)

References
-----------

- `pychem <https://www.astro.unige.ch/~revazy/PyChem/>`_: python module used to generate GEAR tables.
- `Kroupa (2001) <https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..231K/abstract>`_: initial mass function.
- `Poirier (2004) <https://theses.fr/2004STR13003>`_ and `Hausammann (2021) <https://infoscience.epfl.ch/entities/publication/3e6d2e54-a782-440a-86c3-05482e83794d>`_: stellar lifetime fit.
- `Tsujimoto et al. (1995) <https://ui.adsabs.harvard.edu/abs/1995MNRAS.277..945T/abstract>`_ and `Kobayashi et al. (2000) <https://ui.adsabs.harvard.edu/abs/2000ApJ...539...26K/abstract>`_: SNII/SNIa yields.
- `Deng et al. (2024b) <https://arxiv.org/abs/2405.08869>`_: stellar wind model.
