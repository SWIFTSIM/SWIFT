.. GEAR sub-grid model radiation tables
   Darwin Roduit, 24th September 2026

.. _gear_radiation_tables:

Radiation tables
===================

Every radiation channel reads the star's photon output from a ``Data/Radiation`` group in the stellar evolution table (``GEARFeedback:yields_table``, and ``yields_table_first_stars`` for population III stars). A table without this group cannot drive any radiation channel: SWIFT stops at start-up.

Required datasets
--------------------

Every dataset is indexed by the star's mass (a mass-only, "M" table) or by mass and metallicity (a mass x metallicity, "M,Z" table); see "Getting a table" below for which of the published tables is which. SWIFT reads every dataset's own CGS units (a table without a correct ``units`` attribute on a dataset is rejected at load); the units below are what the reader expects on the CGS side, converted to internal units on read:

Whenever any radiation channel is switched on, the group must carry:

* ``Luminosity`` and ``Integrated_Luminosity`` (erg/s and erg/s/Msun of stars formed): bolometric luminosity, used by radiation pressure.
* ``Q_H`` and ``Integrated_Q_H`` (1/s and 1/s/Msun): ionizing photon rate, used by HII photoionization.
* ``DotEExcess`` and ``Integrated_DotEExcess`` (erg/s and erg/s/Msun): excess-photon-energy emission rate above the HI ionization threshold.
* ``MeanPhotonEnergyLW`` and ``Integrated_MeanPhotonEnergyLW`` (both erg): photon-number-weighted mean Lyman-Werner photon energy. Required unconditionally, even for a photoionization-only or radiation-pressure-only run that never touches the ISRF. Unlike the other ``Integrated_`` datasets, this one is not a per-Msun rate: it is an intensive per-photon quantity (:math:`L_\mathrm{LW}/Q_\mathrm{LW}`), so it does not scale with the mass formed.

If the ISRF module is on (``GEARFeedback:with_interstellar_radiation_field: 1``), the group must also carry ``L_PE``/``Integrated_L_PE`` and ``L_LW``/``Integrated_L_LW`` (erg/s and erg/s/Msun), the photoelectric and Lyman-Werner band luminosities.

``Teff`` (photospheric effective temperature, K) is optional and has no ``Integrated_`` counterpart; it is read if present and otherwise simply unavailable.

Getting a table
------------------

The tables served by ``examples/GEAR_ICs_and_SCRIPTS/getChemistryTable.sh`` predate the ``Data/Radiation`` group and carry no radiation data at all: they cannot drive any radiation channel.

``examples/GEAR_ICs_and_SCRIPTS/getRadiationTable.sh [table.hdf5]`` fetches a table that does. It looks, in order, for the file already in the current directory, ``$GEAR_RADIATION_TABLE``, pychem's default output location, and finally a published, checksum-verified copy. Three tables are published:

* ``PopII_parsec_spectral.hdf5``: mass x metallicity spectral table, pychem/PARSEC. Carries every dataset this page describes, including the ISRF's ``L_PE``/``L_LW`` bands. **Use this one for a science run.**
* ``PopIII_parsec_spectral.hdf5``: the same, for population III (first) stars.
* ``radiation_fits_popII.hdf5``: a mass-only, blackbody-fits table with no metallicity axis and no ISRF band data. The shipped HII-region and radiation-pressure examples require this exact file: their star masses were chosen to hit this table's own ``Q_H`` at a specific value (for instance, Starbench's 26.75 Msun reproduces Bisbas et al. (2015)'s :math:`10^{49}` photons/s), so swapping in the spectral table would change those examples' checked numbers. It is not a substitute for the spectral table in a real run.

To generate your own table, use pychem's ``pychem_generate_hdf5_parameters`` with the parameter file matching the dimensionality you need (a piecewise-fits parameter file gives a mass-only table, a spectral one gives a mass x metallicity table).

Checking a table
------------------

``examples/GEAR_ICs_and_SCRIPTS/checkRadiationTable.sh <table.h5> [--with-isrf] [--require-1d]`` verifies a table carries what SWIFT needs before you spend time on the rest of an example's setup. Add ``--with-isrf`` to also check the four ISRF band datasets, and ``--require-1d`` if your own analysis script only understands a mass-only table.

Interpolation
---------------

* ``interpolation_size_mass`` (section ``GEARRadiation``): number of points in the mass-axis interpolation grid built from the table. Must be at least 2. Default: 500, matching the shipped tables' own native grid size; a smaller value costs interpolation error near the source fits' own mass-relation breakpoints.

.. code:: YAML

   GEARRadiation:
     interpolation_size_mass: 500   # Number of points in the mass interpolation of the Data/Radiation table (Default: 500)
