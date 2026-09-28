.. GEAR sub-grid model interstellar radiation field
   Darwin Roduit, 21st September 2026

.. _gear_isrf:

Interstellar radiation field
============================

Young stars produce a local, non-ionizing ultraviolet radiation field that heats the gas through the photoelectric effect on dust grains and dissociates molecular hydrogen. The GEAR interstellar radiation field (ISRF) module follows this field in two bands:

- **PE**, 6 to 11.2 eV, which drives photoelectric heating on dust,
- **LW** (Lyman-Werner), 11.2 to 13.6 eV, which drives :math:`\mathrm{H}_2` photodissociation.

Each star carries a band luminosity interpolated from the radiation datasets of its yields table (``GEARFeedback:yields_table``), as a function of its mass, metallicity and age. The luminosity is deposited on the gas particles inside the star's SPH kernel. Every receiving gas particle then attenuates the field it sees by a dust extinction factor :math:`\exp(-\kappa_\mathrm{eff} \Sigma)`, with the column :math:`\Sigma = \rho \, \ell` built from the particle's own density and a path :math:`\ell` set by the receiver-side extinction mechanism (see below; by default, the star-to-particle separation). The dust opacity scales with the particle's metallicity, so a metal-free gas particle is not shielded.

Injection alone illuminates only the stars' immediate neighbourhoods. Optionally, the deposited field is then transported away from the sources with a hyperbolic moment (M1) scheme under a reduced speed of light, so that the radiation reaches gas beyond the source kernels at a finite, controllable propagation speed instead of instantaneously.

The resulting per-particle field is handed to Grackle every step: the PE band sets the photoelectric heating rate, and the LW band sets the :math:`\mathrm{H}_2` photodissociation rate (Grackle mode 2 or 3 only; see :ref:`gear_radiation`'s Grackle cooling mode paragraph).

Working configurations are shipped in ``examples/SubgridTests/StellarFeedback/ISRF/``; each example directory carries its own README with its configure line, run command and check script.

Configuring and compiling
-------------------------

The module is part of the GEAR feedback model and is coupled to Grackle, so a run needs at least:

.. code:: bash

   ./configure --with-feedback=GEAR --with-chemistry=GEAR_10 --with-stars=GEAR \
               --with-cooling=grackle_2 --with-grackle=$GRACKLE_ROOT

See :ref:`gear_radiation` for the per-band Grackle mode gating shared by every radiation channel (and why ``--with-subgrid=GEAR`` is the wrong shortcut here), and :ref:`gear_grackle_cooling` for the cooling module itself.

**The yields table must carry the ISRF datasets.** This module additionally needs the ``L_PE``, ``L_LW``, ``Integrated_L_PE`` and ``Integrated_L_LW`` datasets in the table's ``Data/Radiation`` group, on top of the datasets every radiation channel needs. See :ref:`gear_radiation_tables` for the full requirement and how to get and check a table.

Switching the module on
-----------------------

A single parameter enables the module:

.. code:: YAML

   GEARFeedback:
     with_interstellar_radiation_field: 1   # Master switch of the ISRF module

This switch turns on both channels and forces the Grackle flags they need (the ISRF field, the dust chemistry, the photoelectric heating and, from ``grackle_2`` upwards, the radiative-transfer rate channel) on internally, so you do not set those yourself.

With ``with_interstellar_radiation_field: 1`` and everything else left at its default, stars illuminate their own kernels, the receiving gas is shielded by its own dust column, and no transport takes place. This is enough for the photoelectric heating channel in a well-resolved interstellar medium, and it is the cheapest configuration.

Propagation
-----------

Set ``GEARFeedback:ISRF_propagation: 1`` to transport the injected field away from the sources, at a reduced speed of light :math:`c_\mathrm{hyp}`. This is what makes the transport affordable: the true speed of light would force a prohibitively small time step.

The radiation update rides the gas particle's existing time step rather than running on a separate clock: for the default ``ISRF_c_hyp_scheme`` (``4``) and every scheme except ``2``, :math:`c_\mathrm{hyp}` is derived from whatever time step the particle would already take, so it never constrains it further. Only scheme ``2`` (a fixed fraction of :math:`c`, set by ``ISRF_c_hyp_fixed_fraction_of_c``) fixes the propagation speed first; the particle's time step is then constrained to keep the receiver-side stability condition satisfied, the same way ``HII_rebuild_time_Myr`` constrains a star's time step (see :ref:`gear_radiation_hii`).

``GEARFeedback:ISRF_c_hyp_scheme`` (0 to 4) selects how the propagation speed and the pairwise transport operators are built on each particle. Use the default (``4``) unless you have a specific reason not to. The one alternative worth considering is ``2``, a fixed fraction of the speed of light, useful when you want a propagation speed that does not vary with resolution or time step. The remaining values (``0``, ``1``, ``3``) are earlier variants kept for comparison and are not recommended for a new run. ``ISRF_c_hyp_fixed_fraction_of_c`` must be positive when ``ISRF_c_hyp_scheme`` is ``2`` and must be zero for every other scheme; SWIFT stops at start-up on either mismatch.

``GEARFeedback:ISRF_c_hyp_margin`` sets the stability-margin coefficient that, in turn, sets the radiation time step: larger values propagate the field faster and cost more steps. Its admissible range is tied to the dissipation coefficients below through a joint stability bound, so raising it requires lowering them. SWIFT checks the bound at start-up: if your margin and dissipation coefficients together violate it, the error message reports the exact admissible number for your own settings, rather than a fixed cutoff you would otherwise have to look up.

.. note::
   ``ISRF_propagation`` is a numerical transport model, not a free physical parameter. Leaving it off is a legitimate choice, but a run that turns it on should keep the propagation parameters at their defaults unless it has a specific reason to change them.

Artificial dissipation
----------------------

The hyperbolic update can produce small negative undershoots of the band energy behind a front. A pairwise artificial dissipation, with five parameters (``ISRF_dissipation_alpha_max``, ``ISRF_dissipation_negativity_threshold``, ``ISRF_dissipation_alpha_floor``, ``ISRF_dissipation_floor_h_over_lambda`` and ``ISRF_dissipation_floor_relaxation_residual``), suppresses them; see the complete parameter list below for what each one does. All five default to values suited to a production run. ``ISRF_dissipation_alpha_max`` and ``ISRF_dissipation_alpha_floor`` enter the same joint stability bound as ``ISRF_c_hyp_margin`` above.

Grackle coupling
----------------

Three ``GrackleCooling`` parameters control how the field is turned into heating and dissociation rates. Their defaults are chosen for a run without the ISRF, so a new ISRF run should look at all three.

``GrackleCooling:H2_self_shielding`` (default ``0``)
  How :math:`\mathrm{H}_2` shields itself from the Lyman-Werner field. ``0`` means no shielding, ``2`` sets the shielding column from the kernel support radius, and ``3`` uses Grackle's own local Jeans length. Mode ``1`` is rejected: its length is read from neighbouring Cartesian grid cells, and SWIFT calls Grackle one particle at a time.

  **The default of 0 is the wrong choice for an ISRF run that tracks** :math:`\mathbf{H}_2`. Unshielded, the local Lyman-Werner rate dissociates molecular gas that a real cloud would keep. Set ``3`` unless you have a reason to prefer ``2``; the shipped ``ISRFCosmology`` and ``ISRFH2Photodissociation`` examples both set ``3``. SWIFT emits a start-up warning if the module is on, :math:`\mathrm{H}_2` is tracked and this is still ``0``.

``GrackleCooling:H2_self_shielding_path`` (default ``kernel_diameter``)
  Used by mode ``2`` only. ``kernel_diameter`` sets the shielding column path to twice the kernel support radius, ``kernel_radius`` to one kernel support radius.

``GrackleCooling:photoelectric_heating_efficiency`` (default ``constant``)
  Which photoelectric efficiency Grackle applies to the PE band: ``constant`` (a fixed efficiency of 0.05, Wolfire et al. 1995 Eq. 1), ``wolfire1995`` (the electron-density-dependent efficiency of the same paper, Eq. 2), or ``density_dependent`` (requires a Grackle build that provides it; SWIFT stops at start-up if the module is on and the Grackle it was built against does not). None of the shipped ISRF examples override this: they run at the default, ``constant``.

Two further parameters are worth checking:

- ``GrackleCooling:local_dust_to_gas_ratio`` sets the dust content that the photoelectric rate is proportional to. A negative value uses Grackle's own default.
- ``GrackleCooling:RT_H2_dissociation_rate_cgs`` applies a constant :math:`\mathrm{H}_2` dissociation rate to all gas, on top of the per-particle rate this module computes. Set it to ``0`` in an ISRF run. SWIFT warns if it is nonzero.

Receiver-side extinction
------------------------

``GEARFeedback:ISRF_extinction_path`` (default ``pair_separation``) selects the mechanism that sets the path length :math:`l` of the dust column each gas particle shields itself with. The column is :math:`\Sigma = \rho_j l`, with the receiver's own density.

``constant_kernel_path``
  :math:`l = R\,\gamma_K h_j`, with :math:`R` set by ``ISRF_extinction_path_in_kernel_radii`` (default ``1.0``). :math:`R = 1` is one kernel support radius, the largest path the geometry admits, since the illuminating star sits inside the receiver's own kernel.

``pair_separation``
  :math:`l = r`, the star-to-particle separation of the pair being injected. It is the only mechanism that varies the attenuation across the kernel instead of applying one flat factor, and it is the default.

``temperature_capped_jeans``
  :math:`l = \min(\lambda_J(\min(T, T_\mathrm{cap})), \gamma_K h_j)`, with :math:`T_\mathrm{cap}` set by ``ISRF_extinction_jeans_temperature_cap_K`` (default ``40`` K).

An unrecognised value is a fatal error, not a fallback. This parameter is independent of ``GrackleCooling:H2_self_shielding_path``, since the two shield different processes.

Complete parameter list
-----------------------

The ISRF section of the ``GEARFeedback`` block, with every parameter at its default. The last three are for tests and diagnostics only and must stay at their defaults in a production run:

.. code:: YAML

   GEARFeedback:
     with_interstellar_radiation_field: 0                    # Master switch of the ISRF module
     ISRF_propagation: 0                                     # Transport the injected field with the hyperbolic scheme
     ISRF_extinction_path: pair_separation                   # Receiver-side dust column path mechanism
     ISRF_extinction_path_in_kernel_radii: 1.0               # Path R in kernel support radii, constant_kernel_path only
     ISRF_extinction_jeans_temperature_cap_K: 40             # Jeans-length temperature cap, temperature_capped_jeans only
     ISRF_c_hyp_scheme: 4                                    # Propagation-speed and operator scheme, 0 to 4
     ISRF_c_hyp_margin: 0.5                                  # Stability-margin coefficient C_hyp
     ISRF_c_hyp_fixed_fraction_of_c: 0                       # Reduced speed of light as a fraction of c, scheme 2 only
     ISRF_dissipation_alpha_max: 0.5                         # Ceiling of the triggered dissipation coefficient
     ISRF_dissipation_negativity_threshold: 0.01             # Relative undershoot at which the trigger saturates
     ISRF_dissipation_alpha_floor: 0.5                       # Dissipation floor under the trigger
     ISRF_dissipation_floor_h_over_lambda: 0.5               # Knee of the floor's roll-off
     ISRF_dissipation_floor_relaxation_residual: 0.40        # Relaxation-residual gate of the floor
     ISRF_c_hyp_pin_for_debugging: 0                         # Debugging only, leave at 0
     ISRF_c_hyp_fixed_fraction_timestep_off_for_debugging: 0 # Debugging only, leave at 0
     ISRF_dissipation_alpha_pin_for_debugging: 0             # Debugging only, leave at 0

and the recommended ``GrackleCooling`` entries for an ISRF run tracking :math:`\mathrm{H}_2` (these are not Grackle's own defaults; see :ref:`gear_grackle_cooling` for the full block and its actual defaults):

.. code:: YAML

   GrackleCooling:
     H2_self_shielding: 3                         # Use 3 (local Jeans length) or 2 in an ISRF run tracking H2
     H2_self_shielding_path: kernel_diameter      # Mode 2 only: kernel_diameter or kernel_radius
     photoelectric_heating_efficiency: constant   # constant, wolfire1995 or density_dependent
     local_dust_to_gas_ratio: -1                  # -1 uses Grackle's own default
     RT_H2_dissociation_rate_cgs: 0               # Keep at 0 in an ISRF run

Initial conditions
------------------

Three optional gas fields let an initial-conditions file seed the radiation field directly, for a test setup that starts with a field already in place rather than with a star that produces one. All three are singular, following SWIFT's convention that distinguishes input from output fields, and all three are for test and diagnostic use, not a normal production IC.

.. list-table::
   :header-rows: 1
   :widths: 30 45 25

   * - Name
     - Description
     - Units
   * - ``PESpecificEnergy``
     - Initial specific PE-band energy
     - [U_L^2 U_T^{-2}]
   * - ``LWSpecificEnergy``
     - Initial specific Lyman-Werner-band energy
     - [U_L^2 U_T^{-2}]
   * - ``LWPhotonSpecificEnergy``
     - Initial Lyman-Werner-band photon-number moment, energy-equivalent at a fixed reference photon energy, not a photon count
     - [U_L^2 U_T^{-2}]

An initial-conditions file without any of them is unaffected: every field starts at zero. ``LWSpecificEnergy`` and ``LWPhotonSpecificEnergy`` are coupled at first init, per particle: on any particle where one of the two is zero and the other is not, SWIFT sets the zero one equal to the other, since a Lyman-Werner field with no accompanying photon moment is by definition at the reference photon energy. A particle with a nonzero ``LWSpecificEnergy`` therefore always carries a photon moment.

Snapshot outputs
------------------

The star and gas fields this channel writes are documented on the :ref:`gear_output_isrf` section of the :ref:`gear_output_fields` page.
