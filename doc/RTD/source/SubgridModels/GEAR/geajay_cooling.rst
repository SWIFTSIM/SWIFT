.. geajay cooling module

.. _geajay_cooling:

Cooling: geajay
~~~~~~~~~~~~~~~

geajay is a library for primordial chemistry and radiative cooling. SWIFT can use it in place of
Grackle (:ref:`gear_grackle_cooling`). The two libraries cannot be linked together. The manual of
geajay is the reference for the physics, the options, the tables and the library interface. It is in
the ``doc/`` directory of the geajay repository. This page says how to select the library in SWIFT
and lists the keys.

Configuring and running
+++++++++++++++++++++++

Select the module at configure time. ``N`` is the species network that SWIFT stores for each
particle: 0 (tabulated cooling), 1 (H and He), 2 (with H\ :sub:`2`) or 3 (with HD):

.. code-block:: bash

   ./configure --with-geajay=<prefix of geajay> --with-cooling=geajay_N

Run SWIFT with ``--hydro --cooling``. SWIFT reads the section ``GeajayCooling`` of the parameter
file. Five keys have no default: ``cloudy_table_file``, ``with_uv_background``, ``redshift``,
``with_metal_cooling`` and ``thermal_time_myr``. A key that the chosen options do not read stops
SWIFT at start-up when its value differs from its default, or when it is set at all and has no
default. The message names the key and the reason.

``primordial_chemistry_network`` may select a network below the one chosen at configure time, not
above it.

The keys
++++++++

The tables list the keys of the section ``GeajayCooling`` of ``examples/parameter_example.yml``, with
the wording of that file. The file lists the defaults and shows how to write each key. A key marked
"optional" may be left out.

Table, background and metals
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Key
     - Meaning
   * - ``cloudy_table_file``
     - Name of the Cloudy table.
   * - ``primordial_chemistry_network``
     - (optional) Primordial chemistry network: ``tabulated``, ``H_He``, ``H_He_H2`` or ``H_He_H2_D``
       (or 0 to 3), not above the network selected at configure time (default: that network).
   * - ``with_uv_background``
     - Enable or not the UV background.
   * - ``redshift``
     - Redshift to use (-1 means time based redshift).
   * - ``with_metal_cooling``
     - Enable or not the metal cooling.
   * - ``thermal_time_myr``
     - Time (in Myr) for adiabatic cooling after a feedback event.
   * - ``metal_cooling_cmb_floor``
     - (optional) Enable/disable the CMB temperature floor of the metal cooling (default: 1).
   * - ``agora_cmb_temperature_floor``
     - (optional) Enable/disable the AGORA redshift-dependent CMB temperature floor
       (T\ :sub:`CMB,0` (1+z), applied by SWIFT) (default: 0).
   * - ``table_density_edge``
     - (optional) Cloudy table lookup beyond its n\ :sub:`H` range: ``extrapolate`` (the line of the end
       cell, as Grackle) or ``edge`` (the value at the edge) (default: ``extrapolate``).
   * - ``maximal_density_nucleon_per_cm3``
     - (optional) Largest density passed to geajay, as the nucleon density rho / m\ :sub:`p` [cm\ :sup:`-3`],
       the value of ``GrackleCooling:maximal_density_Hpcm3``. Zero or negative: no cap (default: -1).
   * - ``hydrogen_mass_fraction``
     - (optional) Hydrogen fraction by mass of the metal-free gas. A negative value takes the library's:
       0.76 in networks 1 to 3. In network 0 it is X of ``/Composition`` of the table: the tables of the
       geajay tools are built at X = 0.7614 (He/H 0.0789). Other tables, for example the SWIFT example
       table ``CloudyData_UVB=HM2012.h5``, carry their own composition, and network 0 takes the X of the
       table (0.715768 for a table without a ``/Composition`` group). Network 0 refuses any other value
       (default: -1).

Solver
------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Key
     - Meaning
   * - ``solver``
     - (optional) Solver of networks 1 to 3: 0 cheap (sub-stepped), 1 accurate (BDF) (default: the
       geajay default).
   * - ``max_substeps``
     - (optional) Maximal number of sub-steps of one solve, or 0 for the geajay default: 10000 x 0.1 /
       ``energy_change_limit``, so 33333 at the limit 0.03, or the BDF default with solver 1
       (default: 0).
   * - ``subcycle_on_failure``
     - (optional) Number of retries with a time-step split in smaller pieces after a failed solve
       (0 to 20, default: 2).
   * - ``energy_change_limit``
     - (optional) Largest fraction of the internal energy that one solver sub-step may change, in
       (0, 1], or 0 for the geajay default: 0.03 with ``cheap_energy_update`` 1, 0.1 with 0 (the value
       Grackle fixes). Smaller values take more sub-steps (default: 0).
   * - ``cheap_energy_update``
     - (optional) Energy update of the cheap solver: 0 explicit (as Grackle) or 1 second order
       (default: the geajay default, 1).
   * - ``accurate_relative_tolerance``
     - (optional) Relative tolerance of the accurate solver, in [1e-12, 1e-2] (default: the geajay
       default, 1e-6).
   * - ``accurate_absolute_tolerance``
     - (optional) Absolute tolerance of the accurate solver, as a fraction of the element mass density,
       in [1e-30, 1e-6] (default: the geajay default, 1e-18).

Self-shielding
--------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Key
     - Meaning
   * - ``uv_self_shielding``
     - (optional) Self shielding: ``none``, ``HI``, ``HI_HeI``, ``HI_HeI_no_HeII`` or
       ``HI_HeI_HeII_column`` (the library's methods, or 0 to 4), or ``density_threshold`` (or -1) for
       the threshold below (default: ``none``).
   * - ``uv_shielding_table_file``
     - Required with ``HI_HeI_HeII_column``. Column file of the UV background. No default. The column
       keys are read only with ``HI_HeI_HeII_column``; setting one to a non-default value otherwise
       stops the run. The method needs a species network (``H_He`` or above) and a UV background:
       ``with_uv_background`` 1, or ``ionising_background_model`` ``uv_table`` under
       ``background_mode`` ``per_band``. The factors come from the column file, built from the
       spectrum of the UV rates file in use; geajay refuses a file of another spectrum or redshift
       grid. A restart takes these keys from the restart files, not from the parameter file, and opens
       the column file again at the same path (a relative path must hold from the directory of the
       restart).
   * - ``uv_shielding_length_method``
     - Required with ``HI_HeI_HeII_column``: ``given_length`` (the proper length path * kernel_gamma *
       h * a) or ``Jeans`` (local Jeans length) for the gas; ``none`` is refused there. This method
       applies to gas particles. For a GEAR sink particle the length is its cut-off radius r_cut * a
       (kernel_gamma * h_sink * a with ``use_fixed_cut_off_radius`` 0), whatever this method, times the
       factor below (default: ``none``).
   * - ``uv_shielding_length_factor``
     - (optional) Factor R\ :sub:`sh` on the length, above 0 (default: 1).
   * - ``uv_shielding_path``
     - (optional) Read only with ``uv_shielding_length_method`` ``given_length``. Path in kernel support
       radii: ``kernel_radius`` or ``kernel_diameter`` (default: ``kernel_radius``).
   * - ``uv_shielding_threshold_nucleon_per_cm3``
     - Required only with ``uv_self_shielding`` ``density_threshold``. Nucleon density rho / m\ :sub:`p`
       [cm\ :sup:`-3`] above which the UV background is switched off.
   * - ``h2_shielding_length_method``
     - (optional) H\ :sub:`2` self-shielding length: ``none``, ``given_length`` (the full path below,
       N\ :sub:`H2` = n\ :sub:`H2` * path) or ``Jeans`` (local Jeans length), or 0, 2, 3
       (default: ``none``). Sobolev (1) is not supported.
   * - ``h2_shielding_path``
     - (optional) Path for ``given_length``, also the path of the dust column of the backgrounds:
       ``kernel_radius`` (kernel_gamma * h) or ``kernel_diameter`` (default: ``kernel_radius``).

Backgrounds
-----------

With ``background_mode`` ``grackle`` the cooling passes the library value of
``interstellar_radiation_field`` (1.7 Habing) for every particle. With ``per_band`` it passes 0, the
local field only, to which the background G\ :sub:`0` is added.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Key
     - Meaning
   * - ``background_mode``
     - (optional) Background rules: ``grackle`` (Grackle's) or ``per_band`` (one switch per band; the
       keys below) (default: ``grackle``).
   * - ``ionising_background_model``
     - (optional) Read only with ``background_mode`` ``per_band``. Ionising band: ``default`` (follows
       ``with_uv_background``), ``none`` or ``uv_table`` (default: ``default``).
   * - ``lw_background_model``
     - (optional) Read only with ``background_mode`` ``per_band``. Lyman-Werner band: ``default``,
       ``none``, ``uv_table``, ``constant_j21``, ``Brown+2025_strong`` or ``Brown+2025_weak``; both Brown
       models need an FG20 rates file (default: ``default``).
   * - ``lw_background_j21``
     - (optional) Constant Lyman-Werner background J21 (Grackle ``LWbackground_intensity``); negative:
       the library default. Above 0 it replaces the table k31 and adds to
       ``rt_h2_dissociation_rate_per_s`` (default: -1).
   * - ``isrf_background_model``
     - (optional) Read only with ``background_mode`` ``per_band``. Background G\ :sub:`0`: ``default``,
       ``none``, ``constant`` or ``uv_table`` (``G0_habing`` of a UV rates file, log-linear in 1 + z;
       with ``lw_background_model`` ``Brown+2025_strong`` its 11.2 to 13.6 eV share follows the Brown
       k31, FG20 files only) (default: ``default``).
   * - ``isrf_background_habing``
     - (optional) Read only with ``background_mode`` ``per_band``. Background G\ :sub:`0` of the model
       ``constant`` [Habing] (default: -1, unset).
   * - ``isrf_from_field``
     - (optional) Read G\ :sub:`0` from the particle field (Grackle ``use_isrf_field``) (default: 0).
   * - ``uv_rates_file``
     - (optional) Read only with ``background_mode`` ``per_band``. HDF5 file whose ``/UVBRates`` is the
       UV background, for example an FG20 file; the metal table always comes from ``cloudy_table_file``.
       Leave the line out (an empty value is invalid) to use ``cloudy_table_file``
       (default: ``cloudy_table_file``).
   * - ``ionising_rates_file``
     - (optional) Read only with ``background_mode`` ``per_band``. HDF5 file of the ionising band (k24 to
       k26, k29, k30, photoheating, cross sections). Needs ``ionising_background_model`` ``uv_table``.
       Leave the line out to use ``uv_rates_file`` (default: ``uv_rates_file``).
   * - ``lw_rates_file``
     - (optional) Read only with ``background_mode`` ``per_band``. HDF5 file of the Lyman-Werner band
       (k27, k28, k31). Needs ``lw_background_model`` ``uv_table`` or a Brown model. Leave the line
       out to use ``uv_rates_file`` (default: ``uv_rates_file``).
   * - ``isrf_rates_file``
     - (optional) Read only with ``background_mode`` ``per_band``. HDF5 file of the background G\ :sub:`0`.
       Needs ``isrf_background_model`` ``uv_table``. Leave the line out to use ``uv_rates_file``
       (default: ``uv_rates_file``).
   * - ``dust_background_attenuation``
     - (optional) Read only with ``background_mode`` ``per_band``. Dust attenuation of the backgrounds:
       ``none`` or ``per_band``; ``per_band`` needs ``h2_shielding_length_method`` ``given_length``
       (default: ``none``).
   * - ``dust_cross_section_lw_cm2``
     - (optional) Read only with ``dust_background_attenuation`` ``per_band``. Dust cross-section of the
       LW band (11.2-13.6 eV) [cm\ :sup:`2`] (default: 1.5e-21).
   * - ``dust_cross_section_pe_cm2``
     - (optional) Read only with ``dust_background_attenuation`` ``per_band``. Dust cross-section of the
       PE band (6-11.2 eV) [cm\ :sup:`2`] (default: 8.6e-22).

The dust cross-sections per H nucleon set the dust attenuation exp(-sigma (D/D\ :sub:`sun`) N\ :sub:`H`)
of the LW and PE backgrounds, with N\ :sub:`H` the column of the H\ :sub:`2` shielding. D is the
dust-to-gas ratio of ``dust_to_gas_model``, the one the heating uses. The cross-sections do not change ``uv_self_shielding``. The defaults are those of Kim et al. 2023, Table 3 (WD01 R\ :sub:`V` =
3.1, per H nucleon, mu\ :sub:`H` = 1.4). D/D\ :sub:`sun` uses the reference dust-to-gas ratio 0.009387
(Grackle's) as Kim's Zd' = 1, which is 0.0081 by mass.

Shielded metals with a column-shielded H and He
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A shielded Cloudy table (for example ``HM12_shielded.h5``) holds metal cooling and heating that
already include the self-shielding of the gas. ``uv_self_shielding`` ``HI_HeI_HeII_column`` (the
library method 4, "mode 4") shields only the H I, He I, He II and H\ :sub:`2` rates. Use both
together with these keys:

.. code-block:: yaml

   GeajayCooling:
     background_mode: per_band
     cloudy_table_file: HM12_shielded.h5
     uv_rates_file: HM12_thin.h5
     uv_self_shielding: HI_HeI_HeII_column
     uv_shielding_table_file: HM12_column_shielding.h5
     uv_shielding_length_method: Jeans
     uv_shielding_length_factor: 0.28209479

With ``per_band``, the metals see the shielded table. The H and He photo-rates come from the thin
table ``uv_rates_file``, and the column shields them. The ``/UVBRates`` of the shielded table and of the
thin table are equal, so mode 4 does not shield the H and He rates twice.

Two refusals apply to the other ways of combining these files:

* With ``background_mode`` ``grackle``, geajay takes the UV rates from ``cloudy_table_file``. A
  shielded table records no spectrum hash, so mode 4 refuses it: "the UV rates file ... records no
  spectrum SHA-256 (/Provenance spectrum_sha256; a shielded or old table): no column file can belong
  to it".
* With ``background_mode`` ``grackle``, ``uv_rates_file`` stops the run: "GeajayCooling:uv_rates_file
  is set, but it is not read: GeajayCooling:background_mode is not per_band".

The recipe is tested for H and He (and H\ :sub:`2`) without dust, in one-zone boxes. It is not tested
with the Lyman-Werner band, the background G\ :sub:`0` or ``dust_background_attenuation``.

Heating and rates from the parameter file
-----------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Key
     - Meaning
   * - ``rt_use_fields``
     - (optional) Arrays of ionization and heating rates are provided (default: 0).
   * - ``rt_heating_rate_cgs``
     - (optional) Heating rate in units of / nHI_cgs.
   * - ``rt_hi_ionisation_rate_per_s``
     - (optional) HI ionisation rate [1/s].
   * - ``rt_hei_ionisation_rate_per_s``
     - (optional) HeI ionisation rate [1/s].
   * - ``rt_heii_ionisation_rate_per_s``
     - (optional) HeII ionisation rate [1/s].
   * - ``rt_h2_dissociation_rate_per_s``
     - (optional) H\ :sub:`2` dissociation rate [1/s].
   * - ``volumetric_heating_rates_cgs``
     - (optional) Volumetric heating rate in cgs [erg/s/cm\ :sup:`3`].
   * - ``specific_heating_rates_cgs``
     - (optional) Specific heating rate in cgs [erg/s/g].

Molecules and dust
------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Key
     - Meaning
   * - ``h2_three_body_fit``
     - (optional) H\ :sub:`2` three-body formation rate by name: ``Abel+2002``, ``Palla+1983``,
       ``Cohen+Westberg1983``, ``Flower+Harris2007``, ``Glover2008`` or ``Forrey2013``; 0 to 5 also
       accepted (Grackle) (default: ``Abel+2002``).
   * - ``h2_cie_cooling``
     - (optional) H\ :sub:`2` collision-induced emission cooling: ``none`` or ``Ripamonti+Abel2004``
       (or 0, 1) (default: ``none``).
   * - ``dust_model``
     - (optional) Dust chemistry: ``none`` or ``basic`` (or 0, 1). ``basic`` turns on the gas-grain heat
       exchange, the dust temperature, H\ :sub:`2` formation on grains, photoelectric heating and the
       recombination cooling on grains. It reads the metal mass fraction of each particle to compute
       the dust-to-gas ratio (default: ``none``).
   * - ``dust_to_gas_model``
     - (optional) Dust-to-gas ratio as a function of the metallicity: ``linear`` (D =
       ``reference_dust_to_gas_ratio`` Z / Z\ :sub:`sun`) or ``remy_ruyer_2014`` (broken power law of
       Remy-Ruyer et al. 2014, linear above its transition). Needs a dust switch
       (``dust_model`` ``basic``, ``h2_formation_on_dust`` ``enabled``, ``dust_recombination_cooling``
       ``enabled`` or a ``photoelectric_heating_model`` above ``none``): without one,
       ``remy_ruyer_2014`` stops the run. The library has a third value, ``from_field`` (2), that takes
       the dust density from a particle field. SWIFT has no such field and refuses it
       (default: ``linear``).
   * - ``reference_dust_to_gas_ratio``
     - (optional) Ratio of dust mass to gas mass in the local Universe, above 0 (-1 for the library default, 0.009387; a value of 0 or below other than -1 stops the run).
   * - ``photoelectric_heating_model``
     - (optional) Photoelectric heating of the dust: ``default`` (-1, ``constant_efficiency`` with
       ``dust_model`` ``basic``, else ``none``), ``none`` (0), ``fixed_rate`` (1, the constant
       ``photoelectric_heating_rate``), ``constant_efficiency`` (2, 0.05 G\ :sub:`0` x 1e-24 erg/s per
       H atom, Wolfire et al. 1995 eq. 1), ``grain_charge_epsilon`` (3, an efficiency set by the grain
       charge, Bakes and Tielens 1994 eq. 43 with the coefficients of Wolfire et al. 2003 eq. 24) or
       ``density_epsilon`` (4, an efficiency that depends on the density only, Smith 2026 eqs. A1 and
       A2); 0 to 4 are also accepted. Every model scales with D / D\ :sub:`ref` and the hydrogen density
       and gives no heating above 2e4 K. A model above ``none`` also works with ``dust_model``
       ``none``, because SWIFT passes the metal mass fraction of every particle (default: ``default``).
   * - ``photoelectric_heating_rate``
     - (optional) Heating rate of the model ``fixed_rate`` per H atom at the reference dust-to-gas
       ratio [erg/s]. Read only with ``fixed_rate``. Zero or above, or -1 for the library value,
       8.5e-26 (default: -1).
   * - ``dust_recombination_cooling``
     - (optional) Recombination cooling on grains: ``default`` (-1, on with ``dust_model`` ``basic``,
       else off), ``disabled`` (0) or ``enabled`` (1). It reads the metal mass fraction of each
       particle, so ``enabled`` also works with ``dust_model`` ``none`` (default: ``default``).
   * - ``h2_formation_on_dust``
     - (optional) H\ :sub:`2` formation on dust grains: ``default`` (-1, on with ``dust_model`` ``basic``
       and a network with H\ :sub:`2`, else off), ``disabled`` (0), ``enabled`` (1). ``enabled`` also works
       with ``dust_model`` ``none`` and with the network ``tabulated``. At the network ``H_He``
       (``primordial_chemistry`` 1) the library refuses it with a start-up error, whatever
       ``dust_model`` is. SWIFT does not check this itself (default: ``default``).
   * - ``h2_photodissociation_heating``
     - (optional) Enable/disable the heating by H\ :sub:`2` photodissociation (default: 0).
   * - ``compton_xray_heating``
     - (optional) Enable/disable the Compton exchange with the X-ray background,
       k\ :sub:`X` (T\ :sub:`x` - T) per electron (Madau and Efstathiou 1999). Unlike Grackle's
       ``Compton_xray_heating``, whose term is 1e4 to 1e6 times too large at T = 1e4 K, z = 0 to 3
       (default: 0).

Initial species
---------------

Each key is (optional): the initial number density ratio of one species to H. The value is ignored if
it is -1 (default: -1). The keys are ``initial_hii_to_h_ratio``, ``initial_hei_to_h_ratio``,
``initial_heii_to_h_ratio``, ``initial_heiii_to_h_ratio``, ``initial_di_to_h_ratio``,
``initial_dii_to_h_ratio``, ``initial_hm_to_h_ratio``, ``initial_h2i_to_h_ratio``,
``initial_h2ii_to_h_ratio`` and ``initial_hdi_to_h_ratio``.

Cosmic rays
-----------

With ``cosmic_rays`` ``reduced`` the fields ``CosmicRayValues`` and ``CosmicRayEnergyLossRates`` join
the snapshot, but ``swift -o`` does not list them.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Key
     - Meaning
   * - ``cosmic_rays``
     - (optional) Cosmic-ray network: ``off`` or ``reduced`` (needs ``primordial_chemistry_network``
       ``H_He`` or above) (default: ``off``).
   * - ``cosmic_ray_input``
     - (optional) Read only with ``cosmic_rays`` ``reduced``. Meaning of the per-particle value:
       ``ionisation_rate`` (s\ :sup:`-1` per H), ``energy_density`` (eV cm\ :sup:`-3`) or ``constant``
       (default: ``ionisation_rate``).
   * - ``cosmic_ray_initial_value``
     - Read only with ``cosmic_rays`` ``reduced`` and a per-particle input (``ionisation_rate`` or
       ``energy_density``), then required. Start value of every particle.
   * - ``cosmic_ray_constant_rate``
     - Read only with ``cosmic_ray_input`` ``constant``, then required. Primary ionisation rate per H
       [s\ :sup:`-1`].
   * - ``cosmic_ray_rate_per_energy_density``
     - Read only with ``cosmic_ray_input`` ``energy_density``, then required. Ionisation rate per
       energy density kappa [s\ :sup:`-1` / (eV cm\ :sup:`-3`)].
   * - ``cosmic_ray_coulomb_loss_coefficient``
     - Read only with ``cosmic_ray_input`` ``energy_density``, then required. Coulomb loss coefficient
       [cm\ :sup:`3` s\ :sup:`-1`].
   * - ``cosmic_ray_hadronic_loss_coefficient``
     - Read only with ``cosmic_ray_input`` ``energy_density``, then required. Hadronic loss
       coefficient [cm\ :sup:`3` s\ :sup:`-1`].
   * - ``cosmic_ray_heating``
     - (optional) Read only with ``cosmic_rays`` ``reduced``. Heat per primary ionisation:
       ``atomic_molecular_fit`` (Gong et al. 2017) or ``fixed_energy``
       (default: ``atomic_molecular_fit``).
   * - ``cosmic_ray_heat_per_ionisation``
     - (optional) Read only with ``cosmic_ray_heating`` ``fixed_energy``. Heat per ionisation [eV]
       (default: 20).

Output
++++++

The snapshots hold ``CoolingTimes`` and ``CoolingRates``. To switch them off, set
``Snapshots:select_output_on: 1`` and ``CoolingTimes_Gas: off``, ``CoolingRates_Gas: off`` under
``Default:`` in the output selection file.
