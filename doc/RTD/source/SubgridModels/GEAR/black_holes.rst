.. GEAR sub-grid model black holes
   Darwin Roduit, 5th August 2026

.. _gear_black_holes:

Black holes and AGN feedback
=============================

GEAR's black hole model is adapted from the EAGLE black hole model
(`Booth & Schaye 2009 <https://ui.adsabs.harvard.edu/abs/2009MNRAS.398...53B/abstract>`_,
`Rosas-Guevara et al. 2015 <https://ui.adsabs.harvard.edu/abs/2015MNRAS.454.1038R/abstract>`_,
`Schaye et al. 2015 <https://ui.adsabs.harvard.edu/abs/2015MNRAS.446..521S/abstract>`_)
to use GEAR's own chemistry, tracers, and pressure floor. It handles Bondi-Hoyle
accretion of gas, thermal AGN feedback, black hole repositioning, and black
hole mergers. The model can be selected with the configuration option
``--with-black-holes=GEAR``, or as part of the ``--with-subgrid=GEAR-BH``
meta-option, which also selects the GEAR chemistry, tracers, pressure floor,
stars, star formation, and feedback modules, plus Grackle cooling. The
GEAR black hole module is located in ``src/black_holes/GEAR/``.

Seeding
-------

Black hole particles are converted from gas particles on the fly through
friends-of-friends (FOF) group finding (``src/fof.c``): any FOF group that
does not already contain a black hole and whose total mass exceeds
``FOF:black_hole_seed_halo_mass_Msun`` has its densest gas particle
converted into a black hole seed, with subgrid mass
``GEARAGN:subgrid_seed_mass_Msun``.

Two separate things determine what counts as a "group" and what counts
towards its mass, both configured by the (model-independent)
``FOF:linking_types`` / ``FOF:attaching_types`` arrays:

- **Group membership** is decided only by particles of a *linking* type
  (recommended: dark matter only — ``linking_types: [0, 1, 0, 0, 0, 0, 0]``
  in gas/DM/.../star/black-hole order). Only these particles are used to
  find neighbours and build the FOF graph, so which halo a group
  corresponds to is a DM-only question.
- **Group mass**, which is compared against
  ``FOF:black_hole_seed_halo_mass_Msun``, sums *every* particle attached to
  the group — the linking-type particles themselves plus any particle of
  an *attaching* type (recommended: gas, stars, black holes —
  ``attaching_types: [1, 0, 0, 0, 1, 1, 0]``) that ends up nearest to one
  of them. So the seeding threshold is a total (DM + gas + stars + BH)
  halo mass, not a DM-only mass, even though DM alone decides which
  particles belong to the halo in the first place.

This requires running with both ``--black-holes`` and ``--fof`` (see the
``--help`` output of the ``swift`` executable for the exact combination of
flags, which also depends on whether feedback is enabled). Black hole
particles (``PartType5``) can also be placed directly in the initial
conditions, in which case seeding is skipped for those particles.

At the time of writing, no example in ``examples/`` configures a
``GEARAGN`` section; the closest starting points are the EAGLE examples
under ``examples/EAGLE_ICs/`` and ``examples/IdealisedCluster/`` (FOF-based
seeding in a cosmological or idealised-cluster context) and
``examples/SubgridTests/BlackHoleSwallowing/`` (black holes placed directly
in the initial conditions, no FOF seeding) — all of which use
``EAGLEAGN`` parameters that map directly onto the ``GEARAGN`` names
documented below.

Accretion
---------

The gas accretion rate is computed from the standard Bondi-Hoyle formula,
using either the density- and sound-speed-weighted properties of all gas
neighbours (default), or a separate contribution from each neighbour
(``GEARAGN:use_multi_phase_bondi``). An optional density-dependent boost
factor following `Booth & Schaye (2009)
<https://ui.adsabs.harvard.edu/abs/2009MNRAS.398...53B/abstract>`_ can be
applied (``GEARAGN:with_boost_factor``), as can the angular-momentum-based
suppression term of `Rosas-Guevara et al. (2015)
<https://ui.adsabs.harvard.edu/abs/2015MNRAS.454.1038R/abstract>`_
(``GEARAGN:with_angmom_limiter``). The accretion rate is capped at a
fraction of the Eddington rate, ``GEARAGN:max_eddington_fraction``.

Unlike EAGLE, GEAR does not support the subgrid-Bondi variant of the
accretion model (which relies on an EAGLE-specific subgrid cooling table to
estimate the unresolved gas density and sound-speed around the black hole);
only the dynamical (resolved) gas properties are used.

Gas is accreted either by swallowing whole gas particles
(``GEARAGN:use_nibbling`` set to 0) or by continuously nibbling a small
fraction of the mass of each gas neighbour (``GEARAGN:use_nibbling`` set to
1, the standard EAGLE-style approach). A fraction of the accreted rest-mass
energy, set by ``GEARAGN:radiative_efficiency``, is assumed to be radiated
away rather than added to the black hole's mass.

**Pressure floor vs. entropy floor.** EAGLE overrides the sound-speed of
gas sitting *on the entropy floor* (EAGLE's effective equation of state for
unresolved star-forming gas) with a fixed value, so that gas artificially
kept warm by the floor does not bias the Bondi rate
(``EAGLEAGN:with_fixed_T_near_EoS``). GEAR has no entropy floor — only the
unrelated, purely numerical ``GEARPressureFloor`` (a Jeans-mass resolution
safeguard, not an ISM effective equation of state) — so this switch has no
GEAR equivalent and was removed rather than kept as a no-op. The pressure
floor could not stand in for it either way: the black hole's Bondi
accretion rate uses each gas neighbour's sound-speed computed directly from
its density and internal energy (``hydro_get_comoving_soundspeed``), not
the pressure-floor-adjusted value cached for the hydro force loop
(``p->force.soundspeed``), so the pressure floor never reaches the
accretion-rate calculation.

AGN feedback
------------

Black holes accumulate accreted energy in a reservoir and release it as
thermal energy into surrounding gas particles once the reservoir exceeds
the energy needed to heat a target number of neighbours
(``GEARAGN:AGN_num_ngb_to_heat``) by a temperature increase
``GEARAGN:AGN_delta_T_K`` (optionally varying with black hole mass and
local gas properties, see ``GEARAGN:use_variable_delta_T``). Only a
fraction ``GEARAGN:coupling_efficiency`` of the radiated energy actually
couples to the gas.

Which gas particles receive the energy is decided by
``GEARAGN:AGN_feedback_model``: ``Isotropic`` draws random directions and
heats the neighbours closest to each direction; ``MinimumDistance`` and
``MinimumDensity`` heat the closest or densest neighbours; ``Random`` heats
a pseudo-random subset. All four models use the same ray-based targeting
machinery (``src/rays.h``) as the sink and stellar feedback modules.

Repositioning and mergers
--------------------------

Black holes below ``GEARAGN:max_reposition_mass`` are moved towards the
lowest-potential particle in their kernel each time step
(``GEARAGN:max_reposition_distance_ratio`` limits how far), either
instantaneously or at a prescribed, mass- and density-dependent speed
(``GEARAGN:set_reposition_speed``). Two black holes merge when they are
close enough and gravitationally bound according to
``GEARAGN:merger_threshold_type``; the mass ratio of the merger is recorded
as minor or major using ``GEARAGN:threshold_minor_merger`` and
``GEARAGN:threshold_major_merger``.

There is no dedicated on/off switch for repositioning: it is gated purely
by ``bp->subgrid_mass > GEARAGN:max_reposition_mass``. Since the subgrid
mass is always strictly positive (seeded at
``GEARAGN:subgrid_seed_mass_Msun`` and only growing through accretion),
setting ``GEARAGN:max_reposition_mass: 0`` disables repositioning for every
black hole, for the lifetime of the run. Conversely, an implausibly large
value (e.g. ``1e20``) makes every black hole eligible, regardless of mass —
the convention used in the EAGLE examples.

Chemistry
---------

Each black hole tracks its metal content as element mass fractions
(``struct chemistry_bpart_data``), consistent with the GEAR sink model.
When a black hole is seeded, its fractions are copied from the parent gas
particle's own (unsmoothed) composition. When it accretes — by swallowing
a whole particle, nibbling part of one, or merging with another black hole
— its fractions are updated as the mass-weighted average of its previous
composition and the accreted material's, using the *smoothed* metal mass
fraction of the accreted gas, as is done for sink particles.

Parameters
----------

All black hole parameters live in the ``GEARAGN`` section of the parameter
file. Most are mandatory (no default is assumed); a few are only read when
a related switch is enabled, as noted below.

.. code:: YAML

   GEARAGN:
     # Initialisation
     subgrid_seed_mass_Msun:              1.5e5    # Black hole subgrid mass at creation time, in solar masses.
     use_subgrid_mass_from_ics:           1        # (Optional, default: 1) Use the subgrid mass from the ICs rather than the dynamical mass.
     with_subgrid_mass_check:             1        # (Optional, default: 1) Error out if the ICs' subgrid mass is not strictly positive. Ignored if use_subgrid_mass_from_ics is 0.

     # Accretion
     use_multi_phase_bondi:               0        # Compute the Bondi rate neighbour-by-neighbour instead of from kernel-averaged gas properties.
     with_angmom_limiter:                 0        # Apply the Rosas-Guevara et al. (2015) viscous suppression term.
     viscous_alpha:                       1e6      # Normalisation of the viscous suppression term. Only read if with_angmom_limiter is 1.
     radiative_efficiency:                0.1      # Fraction of the accreted rest-mass energy that is radiated away (must be <= 1).
     max_eddington_fraction:              1.       # Maximal accretion rate, in units of the Eddington rate.
     eddington_fraction_for_recording:    0.1      # Eddington fraction above which the last high-Eddington time is recorded.
     with_boost_factor:                   0        # Apply the Booth & Schaye (2009) density-dependent accretion boost.
     boost_alpha_only:                    0        # Use a constant boost factor rather than the full density-dependent model. Only read if with_boost_factor is 1.
     boost_alpha:                         1.       # Boost factor (constant, or its lowest value in the density-dependent model). Only read if with_boost_factor is 1.
     boost_beta:                          2.       # Power-law slope of the density-dependent boost. Only read if with_boost_factor is 1 and boost_alpha_only is 0.
     boost_n_h_star_H_p_cm3:              0.1      # Normalisation density of the density-dependent boost, in H/cm^3. Only read if with_boost_factor is 1 and boost_alpha_only is 0.
     use_nibbling:                        1        # Continuously nibble a small mass fraction from each gas neighbour, rather than swallowing whole particles.
     min_gas_mass_for_nibbling_Msun:      10       # Minimum gas particle mass allowed after nibbling, in solar masses. Only read if use_nibbling is 1.

     # Feedback
     AGN_feedback_model:                  Isotropic  # One of Random, Isotropic, MinimumDistance, MinimumDensity.
     AGN_use_deterministic_feedback:      1        # Use a deterministic (rather than probabilistic) number of heating events per time step.
     coupling_efficiency:                 0.15     # Fraction of the radiated energy that couples thermally to the gas.
     AGN_delta_T_K:                       3.16228e8  # (Constant, or initial) AGN heating temperature increase, in Kelvin.
     use_variable_delta_T:                0        # Scale the heating temperature with the black hole mass instead of using a constant value.
     AGN_delta_T_mass_norm:               3.16228e8 # Normalisation temperature for the mass-scaling. Only read if use_variable_delta_T is 1.
     AGN_delta_T_mass_reference:          1e8      # Reference black hole mass for the mass-scaling, in solar masses. Only read if use_variable_delta_T is 1.
     AGN_delta_T_mass_exponent:           0.666667 # Power-law exponent of the mass-scaling. Only read if use_variable_delta_T is 1.
     AGN_with_locally_adaptive_delta_T:   0        # Additionally scale the heating temperature with local gas properties. Only read if use_variable_delta_T is 1.
     AGN_delta_T_crit_factor:             1.       # Buffer factor for the numerical-efficiency temperature. Only read if AGN_with_locally_adaptive_delta_T is 1.
     AGN_delta_T_background_factor:       0.5      # Buffer factor for the background temperature. Only read if AGN_with_locally_adaptive_delta_T is 1.
     AGN_delta_T_max:                     3.16228e9 # Maximal heating temperature increase, in Kelvin. Only read if use_variable_delta_T is 1.
     AGN_delta_T_min:                     3.16228e7 # Minimal heating temperature increase, in Kelvin. Only read if use_variable_delta_T is 1.
     AGN_use_nheat_with_fixed_dT:         0        # Base the energy reservoir threshold on the constant AGN_delta_T_K rather than the actual (possibly variable) dT. Only read if use_variable_delta_T is 1.
     AGN_use_adaptive_energy_reservoir_threshold: 0 # Scale the number of neighbours to heat with the accretion rate instead of using a fixed value.
     AGN_nheat_alpha:                     1.       # Normalisation of the adaptive threshold. Only read if AGN_use_adaptive_energy_reservoir_threshold is 1.
     AGN_nheat_maccr_normalisation:       1.       # Reference accretion rate for the adaptive threshold, in solar masses per year. Only read if AGN_use_adaptive_energy_reservoir_threshold is 1.
     AGN_nheat_limit:                     100.     # Hard upper limit on the adaptive threshold. Only read if AGN_use_adaptive_energy_reservoir_threshold is 1.
     AGN_num_ngb_to_heat:                 1.       # (Constant, or initial) target number of gas neighbours to heat per feedback event.

     # Repositioning
     max_reposition_mass:                 1e8      # Maximal black hole (dynamical) mass for repositioning to be applied, in solar masses.
     max_reposition_distance_ratio:       3.0      # Maximal repositioning distance, in units of the gravitational softening.
     with_reposition_velocity_threshold:  0        # Only reposition onto neighbours below a relative velocity threshold.
     max_reposition_velocity_ratio:       1.0      # Maximal relative velocity, in units of the black hole's ambient sound speed. Only read if with_reposition_velocity_threshold is 1.
     min_reposition_velocity_threshold:   -1.0     # Lower bound (in internal units) on the velocity threshold; ignored if negative. Only read if with_reposition_velocity_threshold is 1.
     set_reposition_speed:                0        # Reposition at a prescribed speed instead of instantaneously.
     reposition_coefficient_upsilon:      0.03     # Normalisation of the repositioning speed, in km/s per solar mass. Only read if set_reposition_speed is 1.
     reposition_reference_mass:           1e8      # Reference black hole mass for the repositioning-speed scaling, in solar masses. Only read if set_reposition_speed is 1.
     reposition_exponent_mass:            2.0      # (Optional, default: 2.0) Power-law exponent of the mass-scaling. Only read if set_reposition_speed is 1.
     reposition_reference_n_H:            0.1      # Reference gas density for the repositioning-speed scaling, in H/cm^3. Only read if set_reposition_speed is 1.
     reposition_exponent_n_H:             1.0      # (Optional, default: 1.0) Power-law exponent of the density-scaling. Only read if set_reposition_speed is 1.
     with_potential_correction:           1        # Exclude the black hole's own contribution when comparing gravitational potentials for repositioning.

     # Merging
     threshold_minor_merger:              0.333    # Mass ratio above which a merger is recorded as minor.
     threshold_major_merger:              0.333    # Mass ratio above which a merger is recorded as major.
     merger_threshold_type:               DynamicalEscapeVelocity # One of CircularVelocity, EscapeVelocity, DynamicalEscapeVelocity.
     merger_max_distance_ratio:           3.0      # Maximal merging distance, in units of the gravitational softening.

     # Time-stepping
     minimum_timestep_Myr:                0.1      # (Optional, default: unlimited) Minimum allowed black hole time-step, in Myr.

The generic (model-independent) black hole neighbour-search parameters
(``BlackHoles:resolution_eta``, ``BlackHoles:h_tolerance``, ...) are
optional and default to the hydro scheme's own values. The FOF seeding
parameters (``FOF:black_hole_seed_halo_mass_Msun``, ...) are documented
separately; see :ref:`Fof_Parameter_Description_label`.
