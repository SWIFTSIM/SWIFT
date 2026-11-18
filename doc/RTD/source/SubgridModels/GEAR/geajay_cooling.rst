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
     - (optional) Primordial chemistry network: ``tabulated``, ``H_He``, ``H_He_H2``, ``H_He_H2_D`` or
       ``H_He_H2_D_ions`` (or 0 to 4), not above the network selected at configure time (default: that
       network). ``H_He_H2_D_ions`` adds D-, HD+ and HeH+ (15 species per particle, snapshot fields
       ``DM``, ``HDII`` and ``HeHII``) and needs ``--with-cooling=geajay_4``.
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
       0.76 in networks 1 to 4. In network 0 it is X of ``/Composition`` of the table: the tables of the
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
     - (optional) Solver of networks 1 to 4: 0 cheap (sub-stepped), 1 accurate (BDF) (default: the
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
   * - ``lw_background_sawtooth_suppression``
     - (optional) 0 or 1. Lowers the Lyman-Werner background of the table before reionisation, as in
       Grackle. Needs ``with_uv_background`` 1 and the table k31 (``lw_background_model``
       ``uv_table``); the library stops the run otherwise (default: 0).
   * - ``uv_background_redshift_on``, ``uv_background_redshift_full_on``,
       ``uv_background_redshift_drop``, ``uv_background_redshift_off``
     - (optional) Redshifts where the UV background switches on, is fully on, starts to drop and
       switches off. Read only when a UV table is loaded: ``with_uv_background`` 1, or a ``uv_table`` model under ``background_mode`` ``per_band``. Unset (-99999): the redshift range of the
       Cloudy table, the last redshift for the first two and the first for the others (default: unset).
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
     - (optional) Volumetric heating rate in cgs [erg/s/cm\ :sup:`3`]. A finite number of 0 or more;
       a negative rate (extra cooling) stops the run. Not together with ``specific_heating_rates_cgs``.
   * - ``specific_heating_rates_cgs``
     - (optional) Specific heating rate in cgs [erg/s/g]. A finite number of 0 or more; a negative rate
       stops the run.

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
     - (required with a dust switch, no default) Dust-to-gas ratio as a function of the metallicity:
       ``linear`` (D = ``reference_dust_to_gas_ratio`` Z / Z\ :sub:`sun`) or ``remy_ruyer_2014``
       (broken power law of Remy-Ruyer et al. 2014, linear above its transition). A dust switch is
       ``dust_model`` ``basic``, ``h2_formation_on_dust`` ``enabled``, ``dust_recombination_cooling``
       ``enabled``, a ``photoelectric_heating_model`` above ``none`` or ``dust_background_attenuation``
       ``per_band``. The run stops if a dust switch is on and the key is absent. With no dust switch
       the key may be absent, and ``remy_ruyer_2014`` stops the run. The library does not count the
       attenuation as a dust switch: with the attenuation as the only switch, only ``linear`` starts.
       The library has a third value, ``from_field`` (2), that takes
       the dust density from a particle field. SWIFT has no such field and refuses it.
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
       with ``dust_model`` ``none``. The networks ``tabulated`` (``primordial_chemistry`` 0) and ``H_He``
       (``primordial_chemistry`` 1) have no H\ :sub:`2`: the library refuses the switch there with a
       start-up error, whatever ``dust_model`` is. SWIFT does not check this itself (default: ``default``).
   * - ``h2_photodissociation_heating``
     - (optional) Enable/disable the heating by H\ :sub:`2` photodissociation (default: 0).
   * - ``compton_xray_heating``
     - (optional) Enable/disable the Compton exchange with the X-ray background,
       k\ :sub:`X` (T\ :sub:`x` - T) per electron (Madau and Efstathiou 1999). Unlike Grackle's
       ``Compton_xray_heating``, whose term is 1e4 to 1e6 times too large at T = 1e4 K, z = 0 to 3
       (default: 0).
   * - ``case_b_recombination``
     - (optional) Recombination case of the H and He recombination rates: ``A`` (0) or ``B`` (1)
       (default: ``B``, which is not the library default ``A``).

The adiabatic index that SWIFT was built with (``--with-adiabatic-index``) goes to geajay as
``adiabatic_index``. At 5/3 the library value stays, so a default build changes no bit.

Initial species
---------------

Each key is (optional): the initial number density ratio of one species to H. The value is ignored if
it is -1 (default: -1). The keys are ``initial_hii_to_h_ratio``, ``initial_hei_to_h_ratio``,
``initial_heii_to_h_ratio``, ``initial_heiii_to_h_ratio``, ``initial_di_to_h_ratio``,
``initial_dii_to_h_ratio``, ``initial_hm_to_h_ratio``, ``initial_h2i_to_h_ratio``,
``initial_h2ii_to_h_ratio``, ``initial_hdi_to_h_ratio``, ``initial_dm_to_h_ratio``,
``initial_hdii_to_h_ratio`` and ``initial_hehii_to_h_ratio``. The last three need the network
``H_He_H2_D_ions``. A species without a ratio starts at the mass fraction 1e-20.

Restart files
-------------

A restart file starts with a record of the geajay layout: the network, the number of species per
particle and the sizes of the cooling structures. SWIFT stops with a message if the record differs
from the running build, or if the file has no record. Restart files written before the record
existed (before the network ``H_He_H2_D_ions``) are refused in every network: start again from a
snapshot or from the initial conditions.

With the network ``H_He_H2_D`` (``primordial_chemistry_network`` 3) geajay holds the total deuterium
at a constant D to H mass ratio and rescales D, D\ :sub:`+` and HD to it at every call. With the
network ``H_He_H2_D_ions`` (4) it rescales D\ :sup:`-` and HD\ :sub:`+` too, because each holds one
D atom (HeH\ :sub:`+` holds none). The key
``deuterium_to_hydrogen_mass_ratio`` (optional) sets that ratio. The value -1 is the library value,
6.8e-5 by mass, which is 3.4e-5 by number. A D atom weighs 2 H atoms, so the initial ratios, which are
number ratios to H, give the mass ratio 2 (``initial_di_to_h_ratio`` + ``initial_dii_to_h_ratio`` +
``initial_hdi_to_h_ratio``), plus ``initial_dm_to_h_ratio`` + ``initial_hdii_to_h_ratio`` inside the
bracket with the network 4. An unset ``initial_di_to_h_ratio`` counts as the library total, so when only
another initial D ratio is set the run stops and asks for ``initial_di_to_h_ratio`` too. The run
stops at start-up if that mass ratio differs from ``deuterium_to_hydrogen_mass_ratio`` by more than a
relative 1e-6. For example, ``initial_di_to_h_ratio: 4e-5`` alone needs
``deuterium_to_hydrogen_mass_ratio: 8e-5``. With a network below 3 a value other than -1 stops the run.

Cosmic rays
-----------

With ``cosmic_rays`` ``reduced`` the fields ``CosmicRayValues`` and ``CosmicRayEnergyLossRates`` join
the snapshot, but ``swift -o`` does not list them. ``CosmicRayValues`` is not written with
``cosmic_ray_source`` ``fire2``.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Key
     - Meaning
   * - ``cosmic_rays``
     - (optional) Cosmic-ray network: ``off`` or ``reduced`` (needs ``primordial_chemistry_network``
       ``H_He`` or above) (default: ``off``).
   * - ``cosmic_ray_source``
     - (optional) Read only with ``cosmic_rays`` ``reduced``. Where the cosmic-ray value of a particle
       comes from: ``constant``, ``per_particle`` or ``fire2`` (default: ``constant``). The three modes
       exclude each other; see below.
   * - ``cosmic_ray_input``
     - (optional) Read only with ``cosmic_rays`` ``reduced``. Meaning of the value with the source
       ``constant``: ``ionisation_rate`` (s\ :sup:`-1` per H), ``energy_density`` (eV cm\ :sup:`-3`) or
       ``constant`` (default: ``ionisation_rate``). The sources ``per_particle`` and ``fire2`` use
       ``energy_density``: the key may be absent or ``energy_density``, any other value stops the run.
   * - ``cosmic_ray_initial_value``
     - Read only with ``cosmic_rays`` ``reduced``, the source ``constant`` and a value input
       (``ionisation_rate`` or ``energy_density``), then required. The value of every particle for the
       whole run, in the unit of the input.
   * - ``cosmic_ray_fire2_energy_density``
     - Read only with the source ``fire2``, then required, no default. The fixed cosmic-ray energy
       density e\ :sub:`cr` [eV cm\ :sup:`-3`], finite and above 0.
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

Source of the cosmic-ray value
++++++++++++++++++++++++++++++

``cosmic_ray_source`` chooses one of three modes.

``constant``
   One value for the whole run, ``cosmic_ray_initial_value`` (or ``cosmic_ray_constant_rate`` with the
   input ``constant``). Nothing is stored in the particles; ``CosmicRayValues`` holds the value of the run.

``per_particle``
   Each particle stores a cosmic-ray energy density [eV cm\ :sup:`-3`, proper]. It is 0 at the start:
   no parameter sets it. The cooling only reads it. It never changes the stored value and does not
   subtract the energy that the cosmic rays lose in the cooling. The module that injects the cosmic
   rays (for example supernova feedback) writes the value and depletes it. The module has no decay and
   no transport of the stored value. The value is written to the snapshots as ``CosmicRayValues``.

``fire2``
   The value is computed at each call from the local hydrogen number density,
   e\ :sub:`cr` n\ :sub:`H` / (0.01 + n\ :sub:`H`), with e\ :sub:`cr` from
   ``cosmic_ray_fire2_energy_density`` and the factor of ``geajay_cosmic_ray_density_factor_fire2()``
   (Hopkins et al. 2018, eq. B.25). It is the energy density passed to the library, so the input is
   ``energy_density``. Here n\ :sub:`H` = X (:math:`\rho` - :math:`\rho_Z`) / m\ :sub:`p` [cm\ :sup:`-3`], with
   X the hydrogen mass fraction of the metal-free gas, :math:`\rho` the physical density and :math:`\rho_Z` the
   metal density. Nothing is stored. The value does not depend on the per-particle value and is not
   multiplied by it.

With ``per_particle`` and ``fire2`` the library needs ``cosmic_ray_rate_per_energy_density``,
``cosmic_ray_coulomb_loss_coefficient`` and ``cosmic_ray_hadronic_loss_coefficient`` (input
``energy_density``). ``cosmic_ray_initial_value`` stops the run with these two modes, and
``cosmic_ray_fire2_energy_density`` stops it with the others. The mode ``constant`` is the default.

Functions for the injecting module
++++++++++++++++++++++++++++++++++

These functions exist only with ``--with-cooling=geajay_N``: a module that calls them needs
``#ifdef COOLING_GEAJAY``. They are declared in ``src/cooling/geajay/cooling.h`` and take the
``struct xpart`` of the particle.

.. code-block:: c

   void cooling_set_cosmic_ray_energy_density(struct xpart *xp, double energy_density);
   double cooling_get_cosmic_ray_energy_density(const struct xpart *xp);
   double cooling_get_cosmic_ray_energy_loss_rate(const struct xpart *xp);

``cooling_set_cosmic_ray_energy_density()`` writes the stored energy density [eV cm\ :sup:`-3`, proper].
The value must be finite, at least 0 and fit in a float. The run stops if the source is not
``per_particle``: the value would never be read. ``cooling_cosmic_ray_source_is_per_particle()`` tells
which source the run has. The cooling reads the value at its next call. The stored value is part of the
restart files but not of the initial conditions: a run started from a snapshot begins at 0.

``cooling_get_cosmic_ray_energy_density()`` returns the stored value.

``cooling_get_cosmic_ray_energy_loss_rate()`` returns the energy that the cosmic rays lost per volume
and time in the last solve of the particle, Coulomb and hadronic losses, [eV cm\ :sup:`-3` s\ :sup:`-1`,
proper], averaged over the step. It is 0 before the first solve, with an input other than
``energy_density`` and after a failed solve, which skips the step. The snapshot field
``CosmicRayEnergyLossRates`` holds the same number, except NaN after a failed solve. The cooling does not
subtract it from the stored density.

A feedback module uses them in its own time-step loop. At an injection it adds the new energy density:

.. code-block:: c

   const double e_old = cooling_get_cosmic_ray_energy_density(xp);
   cooling_set_cosmic_ray_energy_density(xp, e_old + e_injected);

After a cooling step it depletes the density by what was used. The step ``dt_s`` is the one that the
cooling averaged the loss over, in seconds: the argument ``dt`` of ``cooling_cool_part()``, which is the
physical step in a cosmological run (not ``dt_therm``), times
``units_cgs_conversion_factor(us, UNIT_CONV_TIME)``. The setter MUST NOT run while the cooling task of the
same particle runs: the stored value has no lock. The injector can also apply its own decay or transport:

.. code-block:: c

   if (cooling_cosmic_ray_source_is_per_particle()) {
     const double loss = cooling_get_cosmic_ray_energy_loss_rate(xp);
     const double e = cooling_get_cosmic_ray_energy_density(xp);
     cooling_set_cosmic_ray_energy_density(xp, fmax(0., e - loss * dt_s));
   }

Output
++++++

The snapshots hold ``CoolingTimes`` and ``CoolingRates``. To switch them off, set
``Snapshots:select_output_on: 1`` and ``CoolingTimes_Gas: off``, ``CoolingRates_Gas: off`` under
``Default:`` in the output selection file.
