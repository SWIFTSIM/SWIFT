.. AGN spin and jet model
   Filip Husko, 1 April 2022

.. AGN_spin_jet:

Model parameters
----------------

Below we give an example of parameter choices applicable for e.g. a 50 Mpc box. The new parameters are from ``include_jets`` and below. Their descriptions are given next to the parameters.

.. code:: YAML

    SPINJETAGN:
        subgrid_seed_mass_Msun:             1e5        # Black hole subgrid mass at creation time in solar masses.
        use_subgrid_mass_from_ics:          1          # (Optional) Use subgrid masses specified in ICs [1, default], or initialise them to particle masses [0]?
        with_subgrid_mass_check:            1          # (Optional) Verify that initial black hole subgrid masses are positive [1, default]. Only used if use_subgrid_mass_from_ics is 1.
        use_multi_phase_bondi:              0          # Compute Bondi rates per neighbour particle [1] or for the smoothed ambient gas around the black hole [0]?
        use_subgrid_gas_properties:         1          # Use subgrid density [1] or dynamical density [0] to calculate BH accretion rates?
        use_krumholz:                       1          # Use Krumholz et al. (2006) [1] or standard Bondi-Hoyle-Lyttleton formula [0] for black hole accretion rates? Only used if multi_phase_bondi is 0.
        with_krumholz_vorticity:            0          # Include the vorticity term in Krumholz et al. formula? Only used if use_multi_phase_bondi is 0.
        with_angmom_limiter:                0          # Are we applying the Rosas-Guevara (2015) viscous time-scale reduction term?
        viscous_alpha:                      1e6        # Normalisation constant of the viscous time-scale in the accretion reduction term. Only used if with_angmom_limiter is 1.
        with_boost_factor:                  0          # Are we using the model from Booth, Schaye (2009)?
        boost_alpha:                        1.         # Lowest value for the accretion effeciency for the Booth, Schaye 2009 accretion model.
        boost_beta:                         2.         # Slope of the power law for the Booth, Schaye 2009 model, set beta to zero for constant alpha models.
        boost_n_h_star_cm3:                 0.1        # Normalization of the power law for the Booth Schaye 2009 model in cgs (cm^-3).
        eddington_fraction_for_recording:   0.1        # Record the last time BHs reached an Eddington ratio above this threshold.
        use_nibbling:                       1          # Continuously transfer small amounts of mass from all gas neighbours to a black hole [1] or stochastically swallow whole gas particles [0]? 
        min_gas_mass_for_nibbling_Msun:     9e5        # Minimum mass for a gas particle to be nibbled from [M_Sun]. Only used if use_nibbling is 1.
        coupling_efficiency:                0.15       # Fraction of the radiated energy that couples to the gas in feedback events.
        AGN_delta_T_K:                      1e8        # Change in temperature to apply to the gas particle in an AGN feedback event in Kelvin.
        AGN_num_ngb_to_heat:                1.         # Target number of gas neighbours to heat in an AGN feedback event.
        with_potential_correction:          1          # Subtract BH's own contribution to the potential of neighbours when determining repositioning targets.
        max_reposition_mass_Msun:           2e8        # Maximal BH mass considered for BH repositioning in solar masses.
        max_reposition_distance_ratio:      3.0        # Maximal distance a BH can be repositioned, in units of the softening length.
        with_reposition_velocity_threshold: 0          # Should we only reposition to particles that move slowly w.r.t. the black hole?
        max_reposition_velocity_ratio:      0.25       # Maximal velocity offset of a particle to reposition a BH to, in units of the ambient sound speed of the BH. Only meaningful if with_reposition_velocity_ratio is 1.
        min_reposition_velocity_threshold_km_p_s: -1.0 # Minimal value of the velocity threshold for repositioning [km/s], set to < 0 for no effect. Only meaningful if with_reposition_velocity_ratio is 1.
        set_reposition_speed:               0          # Should we reposition black holes with (at most) a prescribed speed towards the potential minimum?
        reposition_coefficient_upsilon:     0.001      # Repositioning speed normalisation [km/s/M_sun]. Only meaningful if set_reposition_speed is 1.
        reposition_exponent_xi:             1.0        # (Optional) Scaling of repositioning velocity with BH subgrid mass (default: 1.0, linear). Only meaningful if set_reposition_speed is 1.
        threshold_major_merger:             0.333      # Mass ratio threshold to consider a BH merger as 'major'
        threshold_minor_merger:             0.1        # Mass ratio threshold to consider a BH merger as 'minor'
        merger_threshold_type:              DynamicalEscapeVelocity   # Type of velocity threshold for BH mergers (CircularVelocity as in EAGLE, EscapeVelocity, or DynamicalEscapeVelocity)
        merger_max_distance_ratio:          3.0        # Maximal distance over which two BHs can merge, in units of the softening length.
        AGN_use_deterministic_feedback:     1          # Deterministic (1) or stochastic (0) AGN feedback model
        AGN_feedback_model:                 Isotropic  # AGN feedback model (Isotropic or MinimumDistance)
        minimum_timestep_yr:                1000.0     # Minimum time-step of black-hole particles
        max_eddington_fraction:             1.         # Maximal allowed accretion rate in units of the Eddington rate.
        include_jets:                       1          # Global switch whether to include jet feedback [1] or not [0].
        turn_off_radiative_feedback:        0          # Global switch whether to turn off radiative (thermal) feedback [1] or not [0]. This should only be used if 'include_jets' is set to 1, since we want feedback in some form or another.
        alpha_acc:                          0.2        # Viscous alpha of the subgrid accretion disks. Likely to be within the 0.1-0.3 range. The main effect is that it sets the transition accretion rate between the thin and thick disk, as dot(m) = 0.2 * alpha^2.
        mdot_crit_ADAF:                     0.01       # The transition normalized accretion rate (Eddington ratio) at which the disc goes from thick (low accretion rates) to thin (high accretion rates). The feedback also changes from kinetic jets to thermal isotropic, respectively.
        seed_spin:                          0.01       # The (randomly-directed) black hole spin assigned to BHs when they are seeded. Should be strictly between 0 and 1.
        AGN_jet_velocity_model:             Constant   # How AGN jet velocities are calculated. If 'Constant', a single value is used. If 'BlackHoleMass', then an empirical relation between halo mass and black hole mass is used to calculate jet velocities. 'HaloMass' is currently not supported. 
        v_jet_km_p_s:                       3160.      # Jet velocity to use if 'AGN_jet_velocity_model' is 'Constant'. Units are km/s.
        v_jet_BH_mass_scaling_reference_mass_Msun: 1e9 # The reference mass used in the relation between halo mass and BH mass used to calculate jet velocities. Only used if 'AGN_jet_velocity_model' is 'BlackHoleMass'.
        v_jet_BH_mass_scaling_slope:        0.5        # The slope of the relation between halo mass and BH mass used to calculate jet velocities. Only used if 'AGN_jet_velocity_model' is 'BlackHoleMass'.
        v_jet_min_km_p_s:                   500        # The minimal jet velocity. This is used if 'AGN_jet_velocity_model' is 'BlackHoleMass', 'MassLoading' or 'Local'.
        v_jet_max_km_p_s:                   1e4        # The maximal jet velocity. This is used if 'AGN_jet_velocity_model' is 'BlackHoleMass', 'MassLoading' or 'Local'.
        opening_angle_in_degrees:           7.5        # The half-opening angle of the jet in degrees. Should use values < 15 unless for tests.
        N_jet:                              2          # Target number of particles to kick as part of a single jet feedback event. Should be a multiple of 2 to ensure approximate momentum conservation (we always kick particles in pairs, one from each 'side' of the BH, relative to the spin vector).
        AGN_jet_feedback_model:             MinimumDistance # Which particles to kick from the black hole smoothing kernels. Should be 'SpinAxis', 'MinimumDistance', 'MaximumDistance' or 'MinimumDensity'
        eps_f_jet:                          1.         # Coupling efficiency for jet feedback. No reason to expect this to be less than 1.
        fix_jet_efficiency:                 0          # Global switch whether to fix jet efficiency to a particular value [1], or use a spin-dependant formula [0].
        jet_efficiency:                     0.1        # The constant jet efficiency used if 'fix_jet_efficiency' is set to 1.
        fix_jet_direction:                  0          # Global switch whether to fix the jet direction to be along the z-axis, instead of along the spin vector.
        accretion_efficiency_mode:          Constant   # How the accretion efficiencies are calculated for the thick accretion disc. If 'Constant', the value of 'accretion_efficiency_thick' will be used. If 'Vari
        able', the accretion efficiency will scale with Eddington ratio.
        accretion_efficiency_thick:         0.01       # The accretion efficiency (suppression factor of the accretion rate) to use in the thick disc (ADAF), to represent the effects of subgrid ADIOS winds that take away most of the mass flowing through the accretion disc.
        accretion_efficiency_slim:          1          # The constant accretion efficiency to use in the slim disc, at super-Eddington rates.
        fix_radiative_efficiency:           0          # Global switch whether to fix the radiative efficiency to a particular value [1], or use a spin-dependant formula [0]. 
        radiative_efficiency:               0.1        # The constant jet efficiency used if 'fix_radiative_efficiency' is set to 1. Otherwise, this value is used to define the Eddington accretion rate.
        TD_region:                          B          # How to treat the subgrid accretion disk if it is thin, according to the Shakura & Sunyaev (1973) model. If set to B, region b will be used. If set to C, region c will be used.
        include_GRMHD_spindown:             1          # Whether to include high jet spindown rates from GRMHD simulations [1], or use an analytical formula that assumes extraction of energy from the rotational mass/energy of the BH.
        delta_ADAF:                         0.2        # Electron heating parameter, which controls the strength of radiative feedback in thick disks. Should be between 0.1 and 0.5. This parameter is only used if turn_off_secondary_feedback is set to 0.
        include_slim_disk:                  1          # Global switch whether to include super-Eddington accretion, modeled as the slim disk. If set to 0, disks will be considered thin even at very large accretion rates.
        use_jets_in_thin_disc:              1          # Whether to use jets alongside radiation in the thin disc at moderate Eddington ratios.
        use_ADIOS_winds:                    1          # Whether to include ADIOS winds in the thick disc as thermal isotropic feedback (same channel as thin disc quasar feedback, but with a different efficiency). 
        slim_disc_wind_factor:              1          # The relative efficiency of slim disc winds at super-Eddington rates. If '1', full winds will be used, while '0' will lead to no winds. Any value in between
        those can also be used. The wind is implemented in the thermal isotropic feedback channel.

Most of these parameters should work well generally, and should not be changed except for tests. We will discuss only some of the more important ones. You can choose whether to have only the thick and thin disk (low and high BH accretion rates, respectively), or you can also include the slim disk at super-Eddington rates with ``include_slim_disk``. You can control what type of feedback you (do not) want with ``include_jets`` and ``turn_off_radiative_feedback``. If you choose to turn off jets, everything will be modeled as a thin disk (regardless of accretion rate), since jets go hand-in-hand with the thick and the slim disk. Similarly, if you turn off radiation, everything will be treated as a thick disk.

If you set ``use_var_v_jet:   0``, you will need to change ``v_jet``, the kicking velocity of particles, depending on what system you are simulating. You should typically choose values at least 10 times larger than the sound speed of the hot gas in your most massive haloes (e.g. 1500 km/s for a MW-type galaxy and 10 000 km/s for a :math:`10^{14}` :math:`\mathrm{M}_\odot` halo). If, on the other hand, you set ``use_var_v_jet:   1``, the launching velocities will vary on their own depending on the typical sound speed (virial velocity) of the hot gas in the haloes. You then need to set ``v_jet_cs_ratio`` to values :math:`\gg1` (10-30 works well) in order to have significant shocking.
