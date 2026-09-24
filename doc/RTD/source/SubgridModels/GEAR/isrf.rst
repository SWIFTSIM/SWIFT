.. GEAR sub-grid model interstellar radiation field
   Darwin Roduit, 21st September 2026

.. _gear_isrf:

Interstellar radiation field
============================

Young stars produce a local, non-ionizing ultraviolet radiation field that heats the gas through the photoelectric effect on dust grains and dissociates molecular hydrogen. The GEAR interstellar radiation field (ISRF) module follows this field in two bands:

- **PE**, 6 to 11.2 eV, which drives photoelectric heating on dust,
- **LW** (Lyman-Werner), 11.2 to 13.6 eV, which drives :math:`\mathrm{H}_2` photodissociation.

Each star carries a band luminosity interpolated from the radiation datasets of its yields table (``GEARFeedback:yields_table``), as a function of its mass, metallicity and age. The luminosity is deposited on the gas particles inside the star's SPH kernel. Every receiving gas particle then attenuates the field it sees by a dust extinction factor :math:`\exp(-\kappa_\mathrm{eff} \Sigma)`, with the column :math:`\Sigma = \rho \, \ell` built from the particle's own density and a path :math:`\ell` set by its smoothing length. The dust opacity scales with the particle's metallicity, so a metal-free gas particle is not shielded.

Injection alone illuminates only the stars' immediate neighbourhoods. Optionally, the deposited field is then transported with a hyperbolic flux-relaxation (M1 moment) scheme under a reduced speed of light, so that the radiation reaches gas beyond the source kernels at a finite, controllable propagation speed instead of instantaneously.

The resulting per-particle field is handed to Grackle every step: the PE band sets the photoelectric heating rate, and the LW band sets the :math:`\mathrm{H}_2` photodissociation rate.

The derivations behind the injection, the extinction and the propagation scheme are given in ``theory/GEAR/Radiation/02_fuv_isrf.tex``. Working configurations are shipped in ``examples/SubgridTests/StellarFeedback/ISRF/``; each example directory carries its own README with its configure line, run command and check script.

Configuring and compiling
-------------------------

The module is part of the GEAR feedback model and is coupled to Grackle, so a run needs at least:

.. code:: bash

   ./configure --with-feedback=GEAR --with-chemistry=GEAR_10 --with-stars=GEAR \
               --with-cooling=grackle_2 --with-grackle=$GRACKLE_ROOT \
               --with-tracers=GEAR

Two configuration choices matter for this module in particular.

**The Grackle mode decides which of the two channels you get.** Photoelectric heating works in every Grackle mode, ``grackle_0`` included. :math:`\mathrm{H}_2` is only followed by the 9-species and 12-species networks, so :math:`\mathrm{H}_2` photodissociation requires ``--with-cooling=grackle_2`` or ``--with-cooling=grackle_3``. A run configured with ``grackle_0`` or ``grackle_1`` silently gets photoelectric heating only, because there is no :math:`\mathrm{H}_2` abundance to dissociate. See :ref:`gear_grackle_cooling` for the cooling module itself.

**The snapshot fields need the GEAR tracers.** The radiation field itself lives on the gas particles whatever the tracers choice, but the snapshot outputs listed below are registered by the GEAR tracers module, so ``--with-tracers=GEAR`` is required to see them. See :ref:`gear_tracers`.

**The yields table must carry the radiation datasets.** The band luminosities are read from the ``Data/Radiation`` group of ``GEARFeedback:yields_table``, and this module needs its ``L_PE``, ``L_LW``, ``Integrated_L_PE`` and ``Integrated_L_LW`` datasets. SWIFT stops at start-up on a table without them. The tables fetched by ``examples/GEAR_ICs_and_SCRIPTS/getChemistryTable.sh`` predate that group, so generate a table with pychem's ``pychem_generate_hdf5_parameters`` instead. The same group must also carry ``MeanPhotonEnergyLW`` and ``Integrated_MeanPhotonEnergyLW``; SWIFT requires these two regardless of whether this module is on, so a photoionization-only or radiation-pressure-only run needs them too. Check any table with ``examples/GEAR_ICs_and_SCRIPTS/checkRadiationTable.sh <table.h5>``, adding ``--with-isrf`` to also check the four band datasets above.

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

Set ``GEARFeedback:ISRF_propagation: 1`` to transport the injected field away from the sources. The scheme evolves the band energy density together with its flux, closed by an M1 moment closure and relaxed towards the local steady state at each step. It propagates at a reduced speed of light :math:`c_\mathrm{hyp}`, which is what makes the scheme affordable: the true speed of light would force a prohibitively small time step.

``GEARFeedback:ISRF_c_hyp_scheme`` selects how :math:`c_{\mathrm{hyp},i}` is set on each particle, and which form the pairwise transport operators take. The five values are:

``0``
  :math:`c_{\mathrm{hyp},i} = \min(C_\mathrm{hyp} h_i / \Delta t_i, c)`, with :math:`\Delta t_i` the particle's own time step.

``1``
  Kernel-local: the same formula, but with the longest time step among the particle and every neighbour in its kernel.

``2``
  Fixed fraction of the speed of light, :math:`c_{\mathrm{hyp},i} = f c`, with :math:`f` set by ``ISRF_c_hyp_fixed_fraction_of_c``.

``3``
  Speed as in ``0``, with every pairwise transport and dissipation operator rebuilt as the per-particle generalisation of the reduced-speed-of-light method.

``4``
  The kernel-local speed of ``1`` feeding the operators of ``3``. This is the default.

The speed and the operator rewrite are two independent axes, but the speed schemes themselves are alternatives and not layers: ``ISRF_c_hyp_fixed_fraction_of_c`` must be positive when ``ISRF_c_hyp_scheme`` is ``2`` and must be zero for every other scheme. SWIFT stops at start-up on either mismatch.

``GEARFeedback:ISRF_c_hyp_margin`` is the coefficient :math:`C_\mathrm{hyp}` in the formulas above, and thereby sets the radiation time step. Larger values propagate the field faster and cost more steps. Its admissible range is tied to the dissipation coefficients below through a joint stability bound, so raising it requires lowering them: at the default dissipation coefficients, the margin cannot exceed about 0.571. SWIFT checks the bound at start-up and stops if it is violated.

.. note::
   ``ISRF_propagation`` is a numerical transport model, not a free physical parameter. Leaving it off is a legitimate choice, but a run that turns it on should keep the propagation parameters at their defaults unless it has a specific reason to change them.

Artificial dissipation
----------------------

The hyperbolic update can produce small negative undershoots of the band energy behind a front. A pairwise artificial dissipation suppresses them. It has two parts: a trigger that responds to a particle undershooting its neighbours, and a floor that is always active in optically thin gas, where the trigger cannot see the positive leading edge of a pulse.

The six parameters are:

``ISRF_dissipation_alpha_max`` (default ``0.5``)
  Ceiling of the triggered dissipation coefficient. ``0`` disables the trigger.

``ISRF_dissipation_negativity_threshold`` (default ``0.01``)
  Relative undershoot below the neighbours' kernel-mean field at which the trigger reaches its ceiling.

``ISRF_dissipation_alpha_floor`` (default ``0.5``)
  Floor applied under the trigger, rolling off once the smoothing length exceeds the local screening length. ``0`` disables the floor.

``ISRF_dissipation_floor_h_over_lambda`` (default ``0.5``)
  Knee of that roll-off, as a ratio of smoothing length to screening length.

``ISRF_dissipation_floor_relaxation_residual`` (default ``0.40``)
  Gates the floor on how far the particle's flux is from the discrete steady state, so that the floor acts on fronts and not on a settled field. ``0`` disables the gate.

``ISRF_dissipation_alpha_pin_for_debugging`` (default ``0``)
  Debugging only. See below.

``ISRF_dissipation_alpha_max`` and ``ISRF_dissipation_alpha_floor`` enter the joint stability bound with ``ISRF_c_hyp_margin`` described above.

Debugging parameters
--------------------

Three parameters exist for tests and diagnostics only and must be left at their defaults in a production run:

- ``ISRF_c_hyp_pin_for_debugging`` (default ``0``) pins every particle's propagation speed to a fixed physical value, bypassing the speed-of-light clamp.
- ``ISRF_c_hyp_fixed_fraction_timestep_off_for_debugging`` (default ``0``) drops the radiation time-step term that keeps the receiver-side stability condition satisfied under scheme ``2``.
- ``ISRF_dissipation_alpha_pin_for_debugging`` (default ``0``) holds the dissipation coefficient of both bands at a fixed value, bypassing the trigger and the floor.

Grackle coupling
----------------

Three ``GrackleCooling`` parameters control how the field is turned into heating and dissociation rates. Their defaults are chosen for a run without the ISRF, so a new ISRF run should look at all three.

``GrackleCooling:H2_self_shielding`` (default ``0``)
  How :math:`\mathrm{H}_2` shields itself from the Lyman-Werner field. ``0`` means no shielding, ``2`` sets the shielding column from the kernel support radius, and ``3`` uses Grackle's own local Jeans length. Mode ``1`` is rejected: its length is read from neighbouring Cartesian grid cells, and SWIFT calls Grackle one particle at a time.

  **The default of 0 is the wrong choice for an ISRF run that tracks** :math:`\mathbf{H}_2`. Unshielded, the local Lyman-Werner rate dissociates molecular gas that a real cloud would keep. Set ``3`` unless you have a reason to prefer ``2``. SWIFT emits a start-up warning if the module is on, :math:`\mathrm{H}_2` is tracked and this is still ``0``.

``GrackleCooling:H2_self_shielding_path`` (default ``kernel_diameter``)
  Used by mode ``2`` only. ``kernel_diameter`` sets the shielding column path to twice the kernel support radius, ``kernel_radius`` to one kernel support radius.

``GrackleCooling:photoelectric_heating_efficiency`` (default ``constant``)
  Which photoelectric efficiency Grackle applies to the PE band.

  - ``constant``: a fixed efficiency of 0.05 (Wolfire et al. 1995, Eq. 1).
  - ``wolfire1995``: the electron-density-dependent efficiency of the same paper (Eq. 2).
  - ``density_dependent``: a density-dependent efficiency, which requires a Grackle build that provides it. SWIFT stops at start-up if the module is on and the Grackle it was built against does not.

``constant`` is a reasonable default for a first run. ``wolfire1995`` is the physically richer choice in gas where the electron fraction is followed, that is, from ``grackle_1`` upwards.

Two further parameters are worth checking:

- ``GrackleCooling:local_dust_to_gas_ratio`` sets the dust content that the photoelectric rate is proportional to. A negative value uses Grackle's own default.
- ``GrackleCooling:RT_H2_dissociation_rate_cgs`` applies a constant :math:`\mathrm{H}_2` dissociation rate to all gas, on top of the per-particle rate this module computes. Set it to ``0`` in an ISRF run. SWIFT warns if it is nonzero.

Receiver-side extinction
------------------------

``GEARFeedback:ISRF_extinction_path`` (default ``pair_separation``) selects the mechanism that sets the path length :math:`l` of the dust column each gas particle shields itself with. The column is :math:`\Sigma = \rho_j l`, with the receiver's own density.

``constant_kernel_path``
  :math:`l = R\,\gamma_K h_j`, with :math:`R` set by ``ISRF_extinction_path_in_kernel_radii`` (default ``1.0``). :math:`R = 1` is one kernel support radius, the largest path the geometry admits, since the illuminating star sits inside the receiver's own kernel. :math:`R = 5/12` is exact in the uniform optically thin limit.

``pair_separation``
  :math:`l = r`, the star-to-particle separation of the pair being injected. Its kernel-weighted mean is exactly the :math:`5/12` above, and it is the only mechanism that varies the attenuation across the kernel instead of applying one flat factor.

``temperature_capped_jeans``
  :math:`l = \min(\lambda_J(\min(T, T_\mathrm{cap})), \gamma_K h_j)`, with :math:`T_\mathrm{cap}` set by ``ISRF_extinction_jeans_temperature_cap_K`` (default ``40`` K).

An unrecognised value is a fatal error, not a fallback. This parameter is independent of ``GrackleCooling:H2_self_shielding_path``, since the two shield different processes.

Two earlier defaults exist, and neither is recoverable without setting the mechanism and the float explicitly. A run archived before this parameter existed recorded no value for it and ran at two kernel support radii: reproduce it with ``constant_kernel_path`` and ``ISRF_extinction_path_in_kernel_radii: 2.0``. A run made while ``constant_kernel_path`` was briefly the default ran at one support radius: reproduce it with the same mechanism and ``1.0``.

Complete parameter list
-----------------------

The ISRF section of the ``GEARFeedback`` block, with every parameter at its default:

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
     ISRF_c_hyp_pin_for_debugging: 0                         # Debugging only
     ISRF_c_hyp_fixed_fraction_timestep_off_for_debugging: 0 # Debugging only
     ISRF_dissipation_alpha_pin_for_debugging: 0             # Debugging only

and the matching ``GrackleCooling`` entries:

.. code:: YAML

   GrackleCooling:
     H2_self_shielding: 3                         # Use 3 (local Jeans length) or 2 in an ISRF run tracking H2
     H2_self_shielding_path: kernel_diameter      # Mode 2 only: kernel_diameter or kernel_radius
     photoelectric_heating_efficiency: constant   # constant, wolfire1995 or density_dependent
     local_dust_to_gas_ratio: -1                  # -1 uses Grackle's own default
     RT_H2_dissociation_rate_cgs: 0               # Keep at 0 in an ISRF run

Initial conditions
------------------

Two optional gas fields let an initial-conditions file seed the radiation field directly, for a setup that starts with a field already in place rather than with a star that produces one. Both are singular, following SWIFT's convention that distinguishes input from output fields.

+---------------------------+---------------------------------------------+---------------------+
| Name                      | Description                                 | Units               |
+===========================+=============================================+=====================+
| ``PESpecificEnergy``      | | Initial specific PE-band energy           | [U_L^2 U_T^{-2}]    |
+---------------------------+---------------------------------------------+---------------------+
| ``LWSpecificEnergy``      | | Initial specific Lyman-Werner-band energy | [U_L^2 U_T^{-2}]    |
+---------------------------+---------------------------------------------+---------------------+

An initial-conditions file without them is unaffected: the field starts at zero.

Snapshot outputs
----------------

These gas fields are written by the GEAR tracers module, so they need ``--with-tracers=GEAR``. All of them are physical quantities with no scale-factor exponent of their own.

+-----------------------------------------------+---------------------------------------------+--------------------------------+
| Name                                          | Description                                 | Units                          |
+===============================================+=============================================+================================+
| ``PESpecificEnergies``                        | | Local specific PE-band field              | [U_L^2 U_T^{-2}]               |
+-----------------------------------------------+---------------------------------------------+--------------------------------+
| ``LWSpecificEnergies``                        | | Local specific Lyman-Werner-band field    | [U_L^2 U_T^{-2}]               |
+-----------------------------------------------+---------------------------------------------+--------------------------------+
| ``PESpecificFluxes``                          | | Tracked specific flux moment, PE band     | [U_L^3 U_T^{-3}]               |
+-----------------------------------------------+---------------------------------------------+--------------------------------+
| ``LWSpecificFluxes``                          | | Tracked specific flux moment, LW band     | [U_L^3 U_T^{-3}]               |
+-----------------------------------------------+---------------------------------------------+--------------------------------+
| ``PESpecificFluxDivergences``                 | | Flux-divergence term of the PE update     | [U_L^2 U_T^{-3}]               |
+-----------------------------------------------+---------------------------------------------+--------------------------------+
| ``LWSpecificFluxDivergences``                 | | Flux-divergence term of the LW update     | [U_L^2 U_T^{-3}]               |
+-----------------------------------------------+---------------------------------------------+--------------------------------+
| ``PEArtificialDissipationCoefficients``       | | Dissipation coefficient, PE band          | [-]                            |
+-----------------------------------------------+---------------------------------------------+--------------------------------+
| ``LWArtificialDissipationCoefficients``       | | Dissipation coefficient, LW band          | [-]                            |
+-----------------------------------------------+---------------------------------------------+--------------------------------+

The flux, flux-divergence and dissipation fields are only meaningful when ``ISRF_propagation`` is on; they stay at zero otherwise.

Six further fields are energy-conservation and undershoot diagnostics. They stay at zero unless SWIFT is configured with ``--enable-debugging-checks``:

+---------------------------------------------+---------------------------------------------------+
| Name                                        | Description                                       |
+=============================================+===================================================+
| ``PEMinimumSpecificEnergies``               | | Most negative PE-band value written since the   |
| ``LWMinimumSpecificEnergies``               | | previous snapshot, 0 if none was negative       |
+---------------------------------------------+---------------------------------------------------+
| ``PECumulativeInjectedSpecificEnergies``    | | Cumulative dose drawn from the source reservoir |
| ``LWCumulativeInjectedSpecificEnergies``    | | since first init                                |
+---------------------------------------------+---------------------------------------------------+
| ``PECumulativeAbsorbedSpecificEnergies``    | | Cumulative energy attributed to dust absorption |
| ``LWCumulativeAbsorbedSpecificEnergies``    | | and to the cosmological redshift term           |
+---------------------------------------------+---------------------------------------------------+

The field plus the absorbed total minus the injected total isolates the transport and dissipation residual, which sums towards zero over the whole particle set.
