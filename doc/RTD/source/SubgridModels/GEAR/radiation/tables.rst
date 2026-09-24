.. GEAR sub-grid model radiation tables
   Darwin Roduit, 24th September 2026

.. _gear_radiation_tables:

Radiation tables
===================

Every radiation channel reads the star's photon output from a ``Data/Radiation`` group in the stellar evolution table (``GEARFeedback:yields_table``, and ``yields_table_first_stars`` for population III stars). A table without this group cannot drive any radiation channel: SWIFT stops at start-up.

Required datasets
--------------------

Whenever any radiation channel is switched on, the group must carry:

* ``Luminosity`` and ``Integrated_Luminosity``: bolometric luminosity, used by radiation pressure.
* ``Q_H`` and ``Integrated_Q_H``: ionizing photon rate, used by HII photoionization.
* ``DotEExcess`` and ``Integrated_DotEExcess``: excess-photon-energy emission rate above the HI ionization threshold.
* ``MeanPhotonEnergyLW`` and ``Integrated_MeanPhotonEnergyLW``: required unconditionally, even for a photoionization-only or radiation-pressure-only run, since SWIFT reports the population's mean Lyman-Werner photon energy at start-up regardless of which channel is active.

If the ISRF module is on (``GEARFeedback:with_interstellar_radiation_field: 1``), the group must also carry ``L_PE``/``Integrated_L_PE`` and ``L_LW``/``Integrated_L_LW``, the photoelectric and Lyman-Werner band luminosities.

``Teff`` (photospheric effective temperature) is optional and has no ``Integrated_`` counterpart; it is read if present and otherwise simply unavailable.

Getting a table
------------------

The tables served by ``examples/GEAR_ICs_and_SCRIPTS/getChemistryTable.sh`` predate the ``Data/Radiation`` group and carry no radiation data at all: they cannot drive any radiation channel.

``examples/GEAR_ICs_and_SCRIPTS/getRadiationTable.sh [table.hdf5]`` fetches a table that does. It looks, in order, for the file already in the current directory, ``$GEAR_RADIATION_TABLE``, pychem's default output location, and finally a published, checksum-verified copy. Three tables are published:

* ``PopII_parsec_spectral.hdf5``: mass x metallicity spectral table, pychem/PARSEC. The default choice for population II stars.
* ``PopIII_parsec_spectral.hdf5``: the same, for population III (first) stars.
* ``radiation_fits_popII.hdf5``: a mass-only, blackbody-fits table. The shipped HII-region and radiation-pressure examples are calibrated against this table's own ``Q_H`` and require this exact file, not the spectral one.

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
