/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 * Copyright (c) 2022 Filip Husko (filip.husko@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_SPIN_JET_BLACK_HOLES_PROPERTIES_H
#define SWIFT_SPIN_JET_BLACK_HOLES_PROPERTIES_H

#include "chemistry.h"
#include "hydro_properties.h"
#include "string.h"

enum AGN_feedback_models {
  AGN_isotropic_model,       /*< Isotropic model of AGN feedback */
  AGN_minimum_distance_model /*< Minimum-distance model of AGN feedback */
};

enum AGN_heating_temperature_models {
  AGN_heating_temperature_constant, /*< Use a constant delta_T */
  AGN_heating_temperature_local,    /*< Variable delta_T */
};

enum AGN_jet_feedback_models {
  AGN_jet_minimum_distance_model, /*< Minimum-distance model of AGN feedback */
  AGN_jet_maximum_distance_model, /*< Maximum-distance model of AGN feedback */
  AGN_jet_spin_axis_model,        /*< Kicking-along-spin
                                      axis model of AGN feedback */
  AGN_jet_minimum_density_model   /*< Minimum-density model of AGN feedback */
};

enum AGN_jet_velocity_models {
  AGN_jet_velocity_constant,     /*< Use a constant jet velocity */
  AGN_jet_velocity_BH_mass,      /*< Scale the jet velocity with BH mass */
  AGN_jet_velocity_mass_loading, /*< Assume constant mass loading */
  AGN_jet_velocity_local, /*< Scale the jet velocity such that particles clear
                              the kernel exactly when a new one is launched,
                              as well as making sure that the kernel never
                              runs out of particles */
  AGN_jet_velocity_halo_mass /*< Scale the jet velocity with halo mass
                                 (NOT WORKING) */
};

enum BH_merger_thresholds {
  BH_mergers_circular_velocity,        /*< v_circ at separation, as in EAGLE */
  BH_mergers_escape_velocity,          /*< v_esc at separation */
  BH_mergers_dynamical_escape_velocity /*< combined v_esc_dyn at separation */
};

enum thin_disc_regions {
  TD_region_B, /*< Region B from Shakura & Sunyaev (1973) */
  TD_region_C  /*< Region C from Shakura & Sunyaev (1973) */
};

enum accretion_efficiency_modes {
  BH_accretion_efficiency_constant, /*< Single number */
  BH_accretion_efficiency_variable  /*< Scaling with Eddington ratio */
};

/**
 * @brief Properties of black holes and AGN feedback in the EAGEL model.
 */
struct black_holes_props {

  /* ----- Basic neighbour search properties ------ */

  /*! Resolution parameter */
  float eta_neighbours;

  /*! Target weightd number of neighbours (for info only)*/
  float target_neighbours;

  /*! Smoothing length tolerance */
  float h_tolerance;

  /*! Tolerance on neighbour number  (for info only)*/
  float delta_neighbours;

  /*! Maximal number of iterations to converge h */
  int max_smoothing_iterations;

  /*! Maximal change of h over one time-step */
  float log_max_h_change;

  /* ----- Initialisation properties  ------ */

  /*! Mass of a BH seed at creation time */
  float subgrid_seed_mass;

  /*! Should we use the subgrid mass specified in ICs? */
  int use_subgrid_mass_from_ics;

  /*! Should we enforce positive subgrid masses initially? */
  int with_subgrid_mass_check;

  /* ----- Properties of the accretion model ------ */

  /*! Calculate Bondi accretion rate based on subgrid properties? */
  int use_subgrid_gas_properties;

  /*! Calculate Bondi accretion rate for individual neighbours? */
  int use_multi_phase_bondi;

  /*! Switch between Bondi [0] or Krumholz [1] accretion rates */
  int use_krumholz;

  /*! In Krumholz mode, should we include the vorticity term? */
  int with_krumholz_vorticity;

  /*! Are we applying the angular-momentum-based multiplicative term from
   * Rosas-Guevara et al. (2015)? */
  int with_angmom_limiter;

  /*! Normalisation of the viscuous angular momentum accretion reduction */
  float alpha_visc;

  /*! Maximal fraction of the Eddington rate allowed. */
  float f_Edd;

  /*! Eddington fraction threshold for recording */
  float f_Edd_recording;

  /*! Switch for the Booth, Schaye 2009 model */
  int with_boost_factor;

  /*! Lowest value of the boost of the Booth, Schaye 2009 model */
  float boost_alpha;

  /*! Power law slope for the boost of the Booth, Schaye 2009 model */
  float boost_beta;

  /*! Normalisation density (internal units) for the boost of the Booth, Schaye
   * 2009 model */
  double boost_n_h_star;

  /*! Switch for nibbling mode */
  int use_nibbling;

  /*! Minimum gas particle mass in nibbling mode */
  float min_gas_mass_for_nibbling;

  /* ---- Properties of the feedback model ------- */

  /*! AGN feedback model: isotropic or minimum distance */
  enum AGN_feedback_models feedback_model;

  /*! Is the AGN feedback model deterministic or stochastic? */
  int AGN_deterministic;

  /*! Feedback coupling efficiency of the black holes. */
  float epsilon_f;

  /*! The type of jet velocity scaling to use. */
  enum AGN_heating_temperature_models AGN_heating_temperature_model;

  /*! Temperature increase induced by AGN feedback (Kelvin),
      in the constant-temperature case */
  float AGN_delta_T_desired;

  /* Numerical factor by which we rescale the variable delta_T formula,
     fiducial value is 1 */
  float delta_T_xi;

  /* The minimum heating temperature to apply in the case of variable feedback,
     expressed in Kelvin */
  float delta_T_min;

  /* The maximum heating temperature to apply in the case of variable feedback,
     expressed in Kelvin */
  float delta_T_max;

  /* Constants used to parametrise the Dalla Vecchia & Schaye (2012)
   condition - Eqn 18. */
  float normalisation_Dalla_Vecchia;
  float ref_ngb_mass_Dalla_Vecchia;
  float ref_density_Dalla_Vecchia;

  /*! Number of gas neighbours to heat in a feedback event */
  float num_ngbs_to_heat;

  /* ---- Properties of the repositioning model --- */

  /*! Maximal mass of BH to reposition */
  float max_reposition_mass;

  /*! Maximal distance to reposition, in units of softening length */
  float max_reposition_distance_ratio;

  /*! Switch to enable a relative velocity limit for particles to which the
   * black holes can reposition */
  int with_reposition_velocity_threshold;

  /*! Maximal velocity offset of particles to which the black hole can
   * reposition, in units of the ambient sound speed of the black hole */
  float max_reposition_velocity_ratio;

  /*! Minimum value of the velocity repositioning threshold */
  float min_reposition_velocity_threshold;

  /*! Switch to enable repositioning at fixed (maximum) speed */
  int set_reposition_speed;

  /*! Normalisation factor for repositioning velocity */
  float reposition_coefficient_upsilon;

  /*! Repositioning velocity scaling with black hole mass */
  float reposition_exponent_xi;

  /*! Correct potential of BH? */
  int correct_bh_potential_for_repositioning;

  /* ---- Properties of the merger model ---------- */

  /*! Mass ratio above which a merger is considered 'minor' */
  float minor_merger_threshold;

  /*! Mass ratio above which a merger is considered 'major' */
  float major_merger_threshold;

  /*! Type of merger threshold */
  enum BH_merger_thresholds merger_threshold_type;

  /*! Maximal distance over which BHs merge, in units of softening length */
  float max_merging_distance_ratio;

  /* ---- Common conversion factors --------------- */

  /*! Conversion factor from temperature to internal energy */
  float temp_to_u_factor;

  /* ---- Black hole time-step properties ---------- */

  /*! -- Minimum allowed time-step of BH in internal units */
  float time_step_min;

  /* ---- Black hole accretion disk parameters ---------- */

  /*! Viscous alpha of the accretion disk, and various factors and powers
      involving it. Here alpha_factor_x refers to alpha raised to the power of
      x, expressed as e.g. 0549 if x is 0.549 (these are the decimal
      expressions for powers of alpha that appear in accretion disk theory).
      alpha_factor_x_inv refers to the inverse of a similar number, or
      equivalently alpha raised to the power of -x. alpha_factor_x_inv_10 is
      alpha * 10 raised to the power of -x, a combination that also often
      appears in the literature.  */
  float alpha_acc;
  float alpha_acc_2;
  float alpha_acc_2_inv;
  float alpha_factor_01;
  float alpha_factor_02;
  float alpha_factor_08;
  float alpha_factor_08_inv;
  float alpha_factor_08_inv_10;
  float alpha_factor_0549;
  float alpha_factor_06222;

  /* BH spin magnitude to assign at seeding */
  float seed_spin;

  /*! Transition accretion rate between thick (ADAF) and thin disk. */
  float mdot_crit_ADAF;

  /*! Gas-to-total pressure ratio of the accretion disk */
  float beta_acc;
  float beta_acc_inv;

  /*! Critical accretion rate that separates two different states in the
      thick disk regime, related to radiation */
  float edd_crit_thick;

  /*! Effective adiabatic index of the accretion disk */
  float gamma_acc;

  /*! Epsilon parameter of the ADAF (thick disk), which appears in
      Narayan & Yi (1994, 1995) model */
  float eps_ADAF;

  /*! Numerical coefficient relating sound speed in ADAF (thick disk) to the
      Keplerian velocity, which appears in Narayan & Yi (1994, 1995) model */
  float cs_0_ADAF;

  /*! Numerical coefficient relating radial velocity in ADAF (thick disk) to
   the Keplerian velocity, which appears in Narayan & Yi (1994, 1995) model */
  float v_0_ADAF;

  /*! Numerical coefficient relating angular velocity in ADAF (thick disk) to
      the Keplerian angular velocity, which appears in Narayan & Yi
      (1994, 1995) model */
  float omega_0_ADAF;

  /*! Aspect ratio of the ADAF (thick disk), which appears in Narayan & Yi
      (1994, 1995) model */
  float h_0_ADAF;
  float h_0_ADAF_2;

  /*! Electron heating parameter of the ADAF (thick disk), which appears in
      Narayan & Yi (1994, 1995) model */
  float delta_ADAF;

  /*! The gamma parameter of the slim disk, which appears in the Wang & Zhou
      (1999) model */
  float gamma_SD;
  float gamma_SD_inv;

  /*! The ratio of vertical to horizontal kinematic viscosity of the thin disk.
      Note that we use the definition xi = nu_2 / nu_1, whereas xi =
      (nu_2 / nu_1) * 2 * alpha^2 is often used, e.g. Fiacconi et al. (2018) */
  float xi_TD;

  /*! Which region of the thin disc (Shakura & Sunyaev) we are assuming the
      subgrid accretion disk is represented as:
       B - region b; gas pressure dominates over radiation pressure,
                     electron scattering dominates the opacity
       C - region c; gas pressure dominates over radiation pressure,
                     free-free absorption dominates the opacity */
  enum thin_disc_regions TD_region;

  /* ---- Jet feedback - related parameters ---------- */

  /*! Global switch for whether to include jets [1] or not [0]. */
  int include_jets;

  /*! Global switch for whether to turn off radiative feedback [1] or not [0].
   */
  int turn_off_radiative_feedback;

  /* Whether we want to include super-Eddington accretion, modeled as the slim
     disk */
  int include_slim_disk;

  /* Whether or not to use jets from the thin disc regime (at moderate
   * Eddington ratios. */
  int use_jets_in_thin_disc;

  /*! Whether to fix the radiative efficiency to some value [1] or not [0]. */
  int fix_radiative_efficiency;

  /*! The radiative efficiency to use if fix_radiative_efficiency is 1. If
       fix_radiative_efficiency is set to 0, this will still be used to
       define the Eddington accretion rate. */
  float radiative_efficiency;

  /*! How many particles to aim to kick as part of jet feedback. */
  int N_jet;

  /*! The type of jet velocity scaling to use. */
  enum AGN_jet_velocity_models AGN_jet_velocity_model;

  /*! Jet velocity if the constant velocity model is used */
  float v_jet;

  /*! Sets the launching velocity of the jet to v_jet_cs_ratio times the
      sound speed of the hot gas in the halo, assuming it is at virial
      temperature. This is used if the launching model is BH_mass or
      halo_mass. */
  float v_jet_cs_ratio;

  /*! The mass loading to use if the launching velocity model is set to use
      a constant mass loading. */
  float v_jet_mass_loading;

  /*! The free numerical parameter to scale the velocity by, if the local or
      sound_speed launching models are used. */
  float v_jet_xi;

  /*! The reference BH mass to use in the case that we employ a BH mass scaling
   * for the jet velocity. */
  float v_jet_BH_mass_scaling_reference_mass;

  /*! The power law slope to use in the case that we employ a BH mass scaling
   * for the jet velocity. */
  float v_jet_BH_mass_scaling_slope;

  /*! The minimal jet velocity to use in the variable-velocity models */
  float v_jet_min;

  /*! The maximal jet velocity to use in the variable-velocity models */
  float v_jet_max;

  /*! The minimum sound speed of hot gas to count when calculating the
   *  smoothed sound speed of gas in the kernel */
  float sound_speed_hot_gas_min;

  /*! The effective (half-)opening angle of the jet. */
  float opening_angle;

  /*! The coupling efficiency for jet feedback. */
  float eps_f_jet;

  /*! Whether to fix the jet efficiency to some value [1] or not [0]. */
  int fix_jet_efficiency;

  /*! The jet efficiency to use if fix_jet_efficiency is 1. */
  float jet_efficiency;

  /* Whether to fix the jet directions to be along the z-axis. */
  int fix_jet_direction;

  /*! The accretion efficiency (suppression of accretion rate) to use in
   *  the thick disc regime (at low Eddington ratios). */
  float accretion_efficiency_thick;

  /*! The accretion efficiency (suppression of accretion rate) to use in
   *  the slim disc regime (at super-Eddington ratios). */
  float accretion_efficiency_slim;

  /*! Expontent to use for scaling of accretion efficiency with transition
   *  radius in the thick disc. */
  float ADIOS_s;

  /* Whether or not we want to use wind feedback in the ADAF/ADIOS regime
     (at low Eddington ratios). */
  int use_ADIOS_winds;

  /* The factor by which we multiply the slim disc wind efficiency - 0 meaning
     no winds, and 1 meaning full winds. */
  float slim_disc_wind_factor;

  /*! The jet launching scheme to use: minimum distance,
      maximum distance, closest to spin axis or minimum density. */
  enum AGN_jet_feedback_models jet_feedback_model;

  /*! The accretion efficiency mode to use: constant or variable
   * (Eddington-ratio dependent) . */
  enum accretion_efficiency_modes accretion_efficiency_mode;
};

/**
 * @brief Initialise the black hole properties from the parameter file.
 *
 * For the basic black holes neighbour finding properties we use the
 * defaults from the hydro scheme if the users did not provide specific
 * values.
 *
 * @param bp The #black_holes_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param hydro_props The already read-in properties of the hydro scheme.
 * @param cosmo The cosmological model.
 */
INLINE static void black_holes_props_init(struct black_holes_props *bp,
                                          const struct phys_const *phys_const,
                                          const struct unit_system *us,
                                          struct swift_params *params,
                                          const struct hydro_props *hydro_props,
                                          const struct cosmology *cosmo) {

  /* Read in the basic neighbour search properties or default to the hydro
     ones if the user did not provide any different values */

  /* Kernel properties */
  bp->eta_neighbours = parser_get_opt_param_float(
      params, "BlackHoles:resolution_eta", hydro_props->eta_neighbours);

  /* Tolerance for the smoothing length Newton-Raphson scheme */
  bp->h_tolerance = parser_get_opt_param_float(params, "BlackHoles:h_tolerance",
                                               hydro_props->h_tolerance);

  /* Get derived properties */
  bp->target_neighbours = pow_dimension(bp->eta_neighbours) * kernel_norm;
  const float delta_eta = bp->eta_neighbours * (1.f + bp->h_tolerance);
  bp->delta_neighbours =
      (pow_dimension(delta_eta) - pow_dimension(bp->eta_neighbours)) *
      kernel_norm;

  /* Number of iterations to converge h */
  bp->max_smoothing_iterations =
      parser_get_opt_param_int(params, "BlackHoles:max_ghost_iterations",
                               hydro_props->max_smoothing_iterations);

  /* Time integration properties */
  const float max_volume_change =
      parser_get_opt_param_float(params, "BlackHoles:max_volume_change", -1);
  if (max_volume_change == -1)
    bp->log_max_h_change = hydro_props->log_max_h_change;
  else
    bp->log_max_h_change = logf(powf(max_volume_change, hydro_dimension_inv));

  /* Initialisation properties  ---------------------------- */

  bp->subgrid_seed_mass =
      parser_get_param_float(params, "SPINJETAGN:subgrid_seed_mass_Msun");

  /* Convert to internal units */
  bp->subgrid_seed_mass *= phys_const->const_solar_mass;

  bp->use_subgrid_mass_from_ics = parser_get_opt_param_int(
      params, "SPINJETAGN:use_subgrid_mass_from_ics", 1);
  if (bp->use_subgrid_mass_from_ics)
    bp->with_subgrid_mass_check = parser_get_opt_param_int(
        params, "SPINJETAGN:with_subgrid_mass_check", 1);

  /* Accretion parameters ---------------------------------- */

  bp->use_subgrid_gas_properties =
      parser_get_param_int(params, "SPINJETAGN:use_subgrid_gas_properties");
  bp->use_multi_phase_bondi =
      parser_get_param_int(params, "SPINJETAGN:use_multi_phase_bondi");
  if (!bp->use_multi_phase_bondi) {
    bp->use_krumholz = parser_get_param_int(params, "SPINJETAGN:use_krumholz");
    bp->with_krumholz_vorticity =
        parser_get_param_int(params, "SPINJETAGN:with_krumholz_vorticity");
  }

  bp->with_angmom_limiter =
      parser_get_param_int(params, "SPINJETAGN:with_angmom_limiter");
  if (bp->with_angmom_limiter)
    bp->alpha_visc = parser_get_param_float(params, "SPINJETAGN:viscous_alpha");

  bp->f_Edd =
      parser_get_param_float(params, "SPINJETAGN:max_eddington_fraction");
  bp->f_Edd_recording = parser_get_param_float(
      params, "SPINJETAGN:eddington_fraction_for_recording");

  /*  Booth Schaye (2009) Parameters */
  bp->with_boost_factor =
      parser_get_param_int(params, "SPINJETAGN:with_boost_factor");

  if (bp->with_boost_factor) {
    bp->boost_alpha = parser_get_param_float(params, "SPINJETAGN:boost_alpha");

    bp->boost_beta = parser_get_param_float(params, "SPINJETAGN:boost_beta");

    /* Load the density in cgs and convert to internal units */
    bp->boost_n_h_star =
        parser_get_param_float(params, "SPINJETAGN:boost_n_h_star_cm3") /
        units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);
  }

  bp->use_nibbling = parser_get_param_int(params, "SPINJETAGN:use_nibbling");
  if (bp->use_nibbling) {
    bp->min_gas_mass_for_nibbling = parser_get_param_float(
        params, "SPINJETAGN:min_gas_mass_for_nibbling_Msun");
    bp->min_gas_mass_for_nibbling *= phys_const->const_solar_mass;
  }

  if ((bp->min_gas_mass_for_nibbling < 1e-5 * bp->subgrid_seed_mass) ||
      (bp->min_gas_mass_for_nibbling > 1e5 * bp->subgrid_seed_mass)) {
    error(
        "The BH seeding mass and minimal gas mass for nibbling differ by more "
        "than 10^5. That is probably indicating a typo in the parameter file.");
  }

  /* Feedback parameters ---------------------------------- */

  char temp[PARSER_MAX_LINE_SIZE];
  parser_get_param_string(params, "SPINJETAGN:AGN_feedback_model", temp);
  if (strcmp(temp, "Isotropic") == 0)
    bp->feedback_model = AGN_isotropic_model;
  else if (strcmp(temp, "MinimumDistance") == 0)
    bp->feedback_model = AGN_minimum_distance_model;
  else
    error(
        "The AGN feedback model must be either MinimumDistance or Isotropic,"
        " not %s",
        temp);

  bp->AGN_deterministic = parser_get_opt_param_int(
      params, "SPINJETAGN:AGN_use_deterministic_feedback", 1);
  bp->epsilon_f =
      parser_get_param_float(params, "SPINJETAGN:coupling_efficiency");

  /* Common conversion factors ----------------------------- */

  /* Calculate temperature to internal energy conversion factor (all internal
   * units) */
  const double k_B = phys_const->const_boltzmann_k;
  const double m_p = phys_const->const_proton_mass;
  const double mu = hydro_props->mu_ionised;
  bp->temp_to_u_factor = k_B / (mu * hydro_gamma_minus_one * m_p);

  /* ---- Black hole time-step properties ------------------ */

  char temp2[PARSER_MAX_LINE_SIZE];
  parser_get_param_string(params, "SPINJETAGN:AGN_delta_T_model", temp2);
  if (strcmp(temp2, "Constant") == 0) {
    bp->AGN_heating_temperature_model = AGN_heating_temperature_constant;

    bp->AGN_delta_T_desired =
        parser_get_param_float(params, "SPINJETAGN:AGN_delta_T_K");
    /* Check that it makes sense. */

    if (bp->AGN_delta_T_desired <= 0.f)
      error("The AGN heating temperature delta T must be > 0 K, not %.5e K.",
            bp->AGN_delta_T_desired);

    bp->AGN_delta_T_desired /=
        units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  } else if (strcmp(temp2, "Local") == 0) {
    bp->AGN_heating_temperature_model = AGN_heating_temperature_local;

    bp->delta_T_xi = parser_get_param_float(params, "SPINJETAGN:delta_T_xi");

    bp->delta_T_min =
        parser_get_param_float(params, "SPINJETAGN:delta_T_min_K");

    /* Check that minimum temperature makes sense */
    if (bp->delta_T_min <= 0.f)
      error(
          "The minimum AGN heating temperature delta T must be > 0 K, "
          "not %.5e K.",
          bp->delta_T_min);

    /* Convert to internal units */
    bp->delta_T_min /= units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

    bp->delta_T_max =
        parser_get_param_float(params, "SPINJETAGN:delta_T_max_K");

    /* Check that minimum temperature makes sense */
    if (bp->delta_T_max <= 0.f)
      error(
          "The maximum AGN heating temperature delta T must be > 0 K, "
          "not %.5e K.",
          bp->delta_T_max);

    /* Convert to internal units */
    bp->delta_T_max /= units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

    float temp_hot_gas_min =
        parser_get_param_float(params, "SPINJETAGN:temperature_hot_gas_min_K");

    temp_hot_gas_min /= units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

    bp->sound_speed_hot_gas_min =
        sqrtf(hydro_gamma * hydro_gamma_minus_one * temp_hot_gas_min *
              bp->temp_to_u_factor);

    /* Define constants used to parametrise the Dalla Vecchia & Schaye (2012)
       condition - Eqn 18. */
    bp->normalisation_Dalla_Vecchia = 1.8e6;
    bp->normalisation_Dalla_Vecchia /=
        units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
    bp->ref_ngb_mass_Dalla_Vecchia = 1e6 * 60. * phys_const->const_solar_mass;

    /* This is nH = 0.1 cm^-3, convert to physical density */
    bp->ref_density_Dalla_Vecchia =
        0.1 * mu * m_p /
        units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);

  } else {
    error(
        "The AGN heating temperature model must be Constant or SoundSpeed,"
        " not %s",
        temp2);
  }

  bp->num_ngbs_to_heat =
      parser_get_param_float(params, "SPINJETAGN:AGN_num_ngb_to_heat");

  /* Reposition parameters --------------------------------- */

  bp->max_reposition_mass =
      parser_get_param_float(params, "SPINJETAGN:max_reposition_mass_Msun");

  /* Convert to internal units */
  bp->max_reposition_mass *= phys_const->const_solar_mass;

  bp->max_reposition_distance_ratio = parser_get_param_float(
      params, "SPINJETAGN:max_reposition_distance_ratio");

  bp->with_reposition_velocity_threshold = parser_get_param_int(
      params, "SPINJETAGN:with_reposition_velocity_threshold");

  if (bp->with_reposition_velocity_threshold) {
    bp->max_reposition_velocity_ratio = parser_get_param_float(
        params, "SPINJETAGN:max_reposition_velocity_ratio");

    /* Prevent nonsensical input */
    if (bp->max_reposition_velocity_ratio <= 0)
      error("max_reposition_velocity_ratio must be positive, not %f.",
            bp->max_reposition_velocity_ratio);

    bp->min_reposition_velocity_threshold = parser_get_param_float(
        params, "SPINJETAGN:min_reposition_velocity_threshold_km_p_s");
    /* Convert from km/s to internal units */
    bp->min_reposition_velocity_threshold *=
        (1e5 / (us->UnitLength_in_cgs / us->UnitTime_in_cgs));
  }

  bp->set_reposition_speed =
      parser_get_param_int(params, "SPINJETAGN:set_reposition_speed");

  if (bp->set_reposition_speed) {
    bp->reposition_coefficient_upsilon = parser_get_param_float(
        params, "SPINJETAGN:reposition_coefficient_upsilon");

    /* Prevent the user from making silly wishes */
    if (bp->reposition_coefficient_upsilon <= 0)
      error(
          "reposition_coefficient_upsilon must be positive, not %f "
          "km/s/M_sun.",
          bp->reposition_coefficient_upsilon);

    /* Convert from km/s to internal units */
    bp->reposition_coefficient_upsilon *=
        (1e5 / (us->UnitLength_in_cgs / us->UnitTime_in_cgs));

    bp->reposition_exponent_xi = parser_get_opt_param_float(
        params, "SPINJETAGN:reposition_exponent_xi", 1.0);
  }

  bp->correct_bh_potential_for_repositioning =
      parser_get_param_int(params, "SPINJETAGN:with_potential_correction");

  /* Merger parameters ------------------------------------- */

  bp->minor_merger_threshold =
      parser_get_param_float(params, "SPINJETAGN:threshold_minor_merger");

  bp->major_merger_threshold =
      parser_get_param_float(params, "SPINJETAGN:threshold_major_merger");

  char temp3[PARSER_MAX_LINE_SIZE];
  parser_get_param_string(params, "SPINJETAGN:merger_threshold_type", temp3);
  if (strcmp(temp3, "CircularVelocity") == 0)
    bp->merger_threshold_type = BH_mergers_circular_velocity;
  else if (strcmp(temp3, "EscapeVelocity") == 0)
    bp->merger_threshold_type = BH_mergers_escape_velocity;
  else if (strcmp(temp3, "DynamicalEscapeVelocity") == 0)
    bp->merger_threshold_type = BH_mergers_dynamical_escape_velocity;
  else
    error(
        "The BH merger model must be either CircularVelocity, EscapeVelocity, "
        "or DynamicalEscapeVelocity, not %s",
        temp3);

  bp->max_merging_distance_ratio =
      parser_get_param_float(params, "SPINJETAGN:merger_max_distance_ratio");

  const double yr_in_cgs = 365.25 * 24. * 3600.;
  bp->time_step_min =
      parser_get_param_float(params, "SPINJETAGN:minimum_timestep_yr") *
      yr_in_cgs / units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /* ---- Black hole accretion disk physics  ------------------ */

  /* The viscosisty parameter of the subgrid accretion disks.*/
  bp->alpha_acc = parser_get_param_float(params, "SPINJETAGN:alpha_acc");

  /* Various factors and powers of alpha, the subgrid accretion disc
     viscosity, that appear in subgrid accretion disk equations. We
     precompute them here. */
  bp->alpha_acc_2 = bp->alpha_acc * bp->alpha_acc;
  bp->alpha_acc_2_inv = 1. / bp->alpha_acc_2;
  bp->alpha_factor_01 = pow(bp->alpha_acc, 0.1);
  bp->alpha_factor_02 = pow(bp->alpha_acc, 0.2);
  bp->alpha_factor_08 = pow(bp->alpha_acc, 0.8);
  bp->alpha_factor_08_inv = 1. / bp->alpha_factor_08;
  bp->alpha_factor_08_inv_10 = pow(bp->alpha_acc * 10., -0.8);
  bp->alpha_factor_0549 = pow(bp->alpha_acc, 0.549);
  bp->alpha_factor_06222 = pow(bp->alpha_acc * 10., 0.6222);

  if ((bp->alpha_acc <= 0.) || (bp->alpha_acc > 1.)) {
    error(
        "The alpha viscosity parameter of accretion disks must be between 0. "
        " and 1., not %f",
        bp->alpha_acc);
  }

  /* The BH seed spin*/
  bp->seed_spin = parser_get_param_float(params, "SPINJETAGN:seed_spin");
  if ((bp->seed_spin <= 0.) || (bp->seed_spin > 1.)) {
    error(
        "The BH seed spin parameter must be strictly between 0 and 1, "
        "not %f",
        bp->seed_spin);
  }

  /* The critical transition accretion rate between the thick and
     thin disk regimes. */
  bp->mdot_crit_ADAF =
      parser_get_param_float(params, "SPINJETAGN:mdot_crit_ADAF");

  /* Calculate the gas-to-total pressure ratio as based on simulations
     (see Yuan & Narayan 2014) */
  bp->beta_acc = 1. / (1. + 2. * bp->alpha_acc);
  bp->beta_acc_inv = 1. / bp->beta_acc;

  /* Calculate the critical accretion rate between two thick disk regimes as
     in Mahadevan (1997). */
  bp->edd_crit_thick = 2. * bp->delta_ADAF * bp->alpha_acc_2 *
                       (1. - bp->beta_acc) * bp->beta_acc_inv;

  /* Calculate the adiabatic index based on how strong the magnetic fields are
    (see Esin 1997) */
  bp->gamma_acc = (8. - 3. * bp->beta_acc) / (6. - 3. * bp->beta_acc);

  /* Calculate numerical factors of the ADAF (thick disk) as in
     Narayan & Yi (1995) */
  bp->eps_ADAF = (1.6667 - bp->gamma_acc) / (bp->gamma_acc - 1.);
  bp->cs_0_ADAF = sqrtf(2. / (5. + 2. * bp->eps_ADAF));
  bp->v_0_ADAF = 3. / (5. + 2. * bp->eps_ADAF);
  bp->omega_0_ADAF = sqrtf(2. * bp->eps_ADAF / (5. + 2. * bp->eps_ADAF));

  /* Instead of using the Narayan & Yi value, we set to 0.3 based on numerous
     GRMHD simulations (see e.g. Narayan et al. 2021) */
  bp->h_0_ADAF = 0.3;
  bp->h_0_ADAF_2 = 0.09;

  bp->delta_ADAF = parser_get_param_float(params, "SPINJETAGN:delta_ADAF");

  if (bp->delta_ADAF <= 0.) {
    error(
        "The delta electron heating parameter of thick accretion disks "
        " must be > 0. not %f",
        bp->delta_ADAF);
  }

  bp->gamma_SD = sqrtf(5.);
  bp->gamma_SD_inv = 1. / bp->gamma_SD;

  /* Formula taken from Lodato et al. (2010) */
  bp->xi_TD = 2. * (1. + 7. * bp->alpha_acc_2) / (4. + bp->alpha_acc_2) /
              bp->alpha_acc_2;

  char temp4[PARSER_MAX_LINE_SIZE];
  parser_get_param_string(params, "SPINJETAGN:TD_region", temp4);
  if (strcmp(temp4, "B") == 0)
    bp->TD_region = TD_region_B;
  else if (strcmp(temp4, "C") == 0)
    bp->TD_region = TD_region_C;
  else
    error("The choice of thin disc region must be B or C, not %s", temp4);

  /* ---- Jet feedback - related parameters ---------- */
  bp->include_jets = parser_get_param_int(params, "SPINJETAGN:include_jets");

  if ((bp->include_jets != 0) && (bp->include_jets != 1)) {
    error("The include_jets parameter must be either 0 or 1, not %d",
          bp->include_jets);
  }

  bp->turn_off_radiative_feedback =
      parser_get_param_int(params, "SPINJETAGN:turn_off_radiative_feedback");

  if ((bp->turn_off_radiative_feedback != 0) &&
      (bp->turn_off_radiative_feedback != 1)) {
    error(
        "The turn_off_radiative_feedback parameter must be either 0 or 1, "
        " not %d",
        bp->turn_off_radiative_feedback);
  }

  if ((bp->turn_off_radiative_feedback) && (bp->include_jets == 0)) {
    error(
        "The turn_off_radiative_feedback parameter and include_jet parameters"
        "cannot at the same time be 1 and 0, respectively. In other words, at"
        "least one of the two feedback modes must be turned on.");
  }

  bp->include_slim_disk =
      parser_get_param_int(params, "SPINJETAGN:include_slim_disk");

  if ((bp->include_slim_disk != 0) && (bp->include_slim_disk != 1)) {
    error(
        "The include_slim_disk parameter must be either 0 or 1, "
        "not %d",
        bp->include_slim_disk);
  }

  bp->use_jets_in_thin_disc =
      parser_get_param_int(params, "SPINJETAGN:use_jets_in_thin_disc");

  if ((bp->use_jets_in_thin_disc != 0) && (bp->use_jets_in_thin_disc != 1)) {
    error(
        "The use_jets_in_thin_disc parameter must be either 0 or 1, "
        "not %d",
        bp->use_jets_in_thin_disc);
  }

  bp->N_jet = parser_get_param_float(params, "SPINJETAGN:N_jet");

  if (bp->N_jet % 2 != 0) {
    error("The N_jet parameter must be divisible by two, not %d", bp->N_jet);
  }

  char temp5[PARSER_MAX_LINE_SIZE];
  parser_get_param_string(params, "SPINJETAGN:AGN_jet_velocity_model", temp5);
  if (strcmp(temp5, "Constant") == 0) {
    bp->AGN_jet_velocity_model = AGN_jet_velocity_constant;

    bp->v_jet = parser_get_param_float(params, "SPINJETAGN:v_jet_km_p_s");

    /* Convert to internal units */
    bp->v_jet *= (1e5 / (us->UnitLength_in_cgs / us->UnitTime_in_cgs));

    if (bp->v_jet <= 0.)
      error("The v_jet parameter must be > 0., not %f", bp->v_jet);

  } else if (strcmp(temp5, "BlackHoleMass") == 0) {
    bp->AGN_jet_velocity_model = AGN_jet_velocity_BH_mass;

    bp->v_jet_BH_mass_scaling_reference_mass = parser_get_param_float(
        params, "SPINJETAGN:v_jet_BH_mass_scaling_reference_mass_Msun");
    bp->v_jet_BH_mass_scaling_reference_mass *= phys_const->const_solar_mass;
    bp->v_jet_BH_mass_scaling_slope = parser_get_param_float(
        params, "SPINJETAGN:v_jet_BH_mass_scaling_slope");

    bp->v_jet_min =
        parser_get_param_float(params, "SPINJETAGN:v_jet_min_km_p_s");
    bp->v_jet_min *= (1e5 / (us->UnitLength_in_cgs / us->UnitTime_in_cgs));

    bp->v_jet_max =
        parser_get_param_float(params, "SPINJETAGN:v_jet_max_km_p_s");
    bp->v_jet_max *= (1e5 / (us->UnitLength_in_cgs / us->UnitTime_in_cgs));

  } else if (strcmp(temp5, "MassLoading") == 0) {
    bp->AGN_jet_velocity_model = AGN_jet_velocity_mass_loading;

    bp->v_jet_mass_loading =
        parser_get_param_float(params, "SPINJETAGN:v_jet_mass_loading");

    bp->v_jet_min =
        parser_get_param_float(params, "SPINJETAGN:v_jet_min_km_p_s");
    bp->v_jet_min *= (1e5 / (us->UnitLength_in_cgs / us->UnitTime_in_cgs));

    bp->v_jet_max =
        parser_get_param_float(params, "SPINJETAGN:v_jet_max_km_p_s");
    bp->v_jet_max *= (1e5 / (us->UnitLength_in_cgs / us->UnitTime_in_cgs));

  } else if (strcmp(temp5, "Local") == 0) {
    bp->AGN_jet_velocity_model = AGN_jet_velocity_local;

    bp->v_jet_xi = parser_get_param_float(params, "SPINJETAGN:v_jet_xi");

    bp->v_jet_min =
        parser_get_param_float(params, "SPINJETAGN:v_jet_min_km_p_s");
    bp->v_jet_min *= (1e5 / (us->UnitLength_in_cgs / us->UnitTime_in_cgs));

    bp->v_jet_max =
        parser_get_param_float(params, "SPINJETAGN:v_jet_max_km_p_s");
    bp->v_jet_max *= (1e5 / (us->UnitLength_in_cgs / us->UnitTime_in_cgs));

    float temp_hot_gas_min =
        parser_get_param_float(params, "SPINJETAGN:temperature_hot_gas_min_K");

    temp_hot_gas_min /= units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

    bp->sound_speed_hot_gas_min =
        sqrtf(hydro_gamma * hydro_gamma_minus_one * temp_hot_gas_min *
              bp->temp_to_u_factor);

  } else if (strcmp(temp5, "HaloMass") == 0) {
    error(
        "The scaling of jet velocities with halo mass is currently not "
        "supported.");
  } else {
    error(
        "The AGN jet velocity model must be Constant, MassLoading, "
        "BlackHoleMass, Local, SoundSpeed or HaloMass, not %s",
        temp5);
  }

  bp->opening_angle =
      parser_get_param_float(params, "SPINJETAGN:opening_angle_in_degrees");
  bp->opening_angle = bp->opening_angle * M_PI / 180.;

  bp->eps_f_jet = parser_get_param_float(params, "SPINJETAGN:eps_f_jet");

  if ((bp->eps_f_jet <= 0.) || (bp->eps_f_jet > 1.)) {
    error(
        "The eps_f_jet corresponding to the jet coupling efficiency "
        "must be between 0. and 1., not %f",
        bp->eps_f_jet);
  }

  bp->fix_jet_efficiency =
      parser_get_param_int(params, "SPINJETAGN:fix_jet_efficiency");

  if ((bp->fix_jet_efficiency != 0) && (bp->fix_jet_efficiency != 1)) {
    error(
        "The fix_jet_efficiency parameter must be either 0 or 1, "
        "not %d",
        bp->fix_jet_efficiency);
  }

  bp->jet_efficiency =
      parser_get_param_float(params, "SPINJETAGN:jet_efficiency");

  if (bp->jet_efficiency <= 0.) {
    error("The constant jet efficiency parameter must be >0., not %f",
          bp->jet_efficiency);
  }

  bp->fix_jet_direction =
      parser_get_param_int(params, "SPINJETAGN:fix_jet_direction");

  if ((bp->fix_jet_direction != 0) && (bp->fix_jet_direction != 1)) {
    error(
        "The fix_jet_direction parameter must be either 0 or 1, "
        "not %d",
        bp->fix_jet_direction);
  }

  bp->fix_radiative_efficiency =
      parser_get_param_int(params, "SPINJETAGN:fix_radiative_efficiency");

  if ((bp->fix_radiative_efficiency != 0) &&
      (bp->fix_radiative_efficiency != 1)) {
    error(
        "The fix_radiative_efficiency parameter must be either 0 or 1, "
        "not %d",
        bp->fix_radiative_efficiency);
  }

  bp->radiative_efficiency =
      parser_get_param_float(params, "SPINJETAGN:radiative_efficiency");

  if ((bp->radiative_efficiency <= 0.) || (bp->radiative_efficiency > 1.)) {
    error(
        "The radiative_efficiency corresponding to the radiative efficiency "
        "must be between 0. and 1., not %f",
        bp->radiative_efficiency);
  }

  bp->use_ADIOS_winds =
      parser_get_param_int(params, "SPINJETAGN:use_ADIOS_winds");

  if ((bp->use_ADIOS_winds != 0) && (bp->use_ADIOS_winds != 1)) {
    error(
        "The use_ADIOS_winds parameter must be either 0 or 1, "
        "not %d",
        bp->use_ADIOS_winds);
  }

  bp->slim_disc_wind_factor =
      parser_get_param_float(params, "SPINJETAGN:slim_disc_wind_factor");

  if ((bp->slim_disc_wind_factor < 0) || (bp->slim_disc_wind_factor > 1)) {
    error(
        "The slim_disc_wind_factor parameter must be between 0 and 1, "
        "(inclusive), not %f",
        bp->slim_disc_wind_factor);
  }

  char temp6[60];
  parser_get_param_string(params, "SPINJETAGN:AGN_jet_feedback_model", temp6);
  if (strcmp(temp6, "MinimumDistance") == 0)
    bp->jet_feedback_model = AGN_jet_minimum_distance_model;
  else if (strcmp(temp6, "MaximumDistance") == 0)
    bp->jet_feedback_model = AGN_jet_maximum_distance_model;
  else if (strcmp(temp6, "SpinAxis") == 0)
    bp->jet_feedback_model = AGN_jet_spin_axis_model;
  else if (strcmp(temp6, "MinimumDensity") == 0)
    bp->jet_feedback_model = AGN_jet_minimum_density_model;
  else
    error(
        "The AGN feedback model must be MinimumDistance, MaximumDistance, "
        "SpinAxis or MinimumDensity, not %s",
        temp6);

  bp->accretion_efficiency_slim =
      parser_get_param_float(params, "SPINJETAGN:accretion_efficiency_slim");

  if ((bp->accretion_efficiency_slim < 0) ||
      (bp->accretion_efficiency_slim > 1)) {
    error(
        "The accretion_efficiency_slim parameter must be between 0 and 1, "
        "(inclusive), not %f",
        bp->accretion_efficiency_slim);
  }

  char temp7[60];
  parser_get_param_string(params, "SPINJETAGN:accretion_efficiency_mode",
                          temp7);
  if (strcmp(temp7, "Constant") == 0) {
    bp->accretion_efficiency_mode = BH_accretion_efficiency_constant;
    bp->accretion_efficiency_thick =
        parser_get_param_float(params, "SPINJETAGN:accretion_efficiency_thick");

    if ((bp->accretion_efficiency_thick < 0) ||
        (bp->accretion_efficiency_thick > 1)) {
      error(
          "The accretion_efficiency_thick parameter must be between 0 and 1, "
          "(inclusive), not %f",
          bp->accretion_efficiency_thick);
    }
  } else if (strcmp(temp7, "Variable") == 0) {
    bp->accretion_efficiency_mode = BH_accretion_efficiency_variable;
    bp->ADIOS_s = parser_get_param_float(params, "SPINJETAGN:ADIOS_s");

    if ((bp->ADIOS_s < 0) || (bp->ADIOS_s > 1)) {
      error(
          "The ADIOS_s parameter must be between 0 and 1, "
          "(inclusive), not %f",
          bp->ADIOS_s);
    }
  } else {
    error("The accretion efficiency model must be Constant or Variable, not %s",
          temp7);
  }
}

/**
 * @brief Write a black_holes_props struct to the given FILE as a stream of
 * bytes.
 *
 * @param props the black hole properties struct
 * @param stream the file stream
 */
INLINE static void black_holes_struct_dump(
    const struct black_holes_props *props, FILE *stream) {
  restart_write_blocks((void *)props, sizeof(struct black_holes_props), 1,
                       stream, "black_holes props", "black holes props");
}

/**
 * @brief Restore a black_holes_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param props the black hole properties struct
 * @param stream the file stream
 */
INLINE static void black_holes_struct_restore(
    const struct black_holes_props *props, FILE *stream) {
  restart_read_blocks((void *)props, sizeof(struct black_holes_props), 1,
                      stream, NULL, "black holes props");
}

#endif /* SWIFT_SPIN_JET_BLACK_HOLES_PROPERTIES_H */
