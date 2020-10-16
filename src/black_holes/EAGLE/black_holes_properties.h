/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_EAGLE_BLACK_HOLES_PROPERTIES_H
#define SWIFT_EAGLE_BLACK_HOLES_PROPERTIES_H

#include "chemistry.h"
#include "hydro_properties.h"

#include <string.h>

enum AGN_feedback_models {
  AGN_random_ngb_model,       /*< Random neighbour model for AGN feedback */
  AGN_isotropic_model,        /*< Isotropic model of AGN feedback */
  AGN_minimum_distance_model, /*< Minimum-distance model of AGN feedback */
  AGN_minimum_density_model   /*< Minimum-density model of AGN feedback */
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

  /*! Calculate Bondi accretion rate for individual neighbours? */
  int use_multi_phase_bondi;

  /*! Are we using the subgrid gas properties in the Bondi model? */
  int use_subgrid_bondi;

  /*! Are we applying the angular-momentum-based multiplicative term from
   * Rosas-Guevara et al. (2015)? */
  int with_angmom_limiter;

  /*! Normalisation of the viscuous angular momentum accretion reduction */
  float alpha_visc;

  /*! Radiative efficiency of the black holes. */
  float epsilon_r;

  /*! Maximal fraction of the Eddington rate allowed. */
  float f_Edd;

  /*! Eddington fraction threshold for recording */
  float f_Edd_recording;

  /*! Switch for the Booth & Schaye 2009 model */
  int with_boost_factor;

  /*! Use constant-alpha version of Booth & Schaye (2009) model? */
  int boost_alpha_only;

  /*! Lowest value of the boost of the Booth & Schaye 2009 model */
  float boost_alpha;

  /*! Power law slope for the boost of the Booth & Schaye 2009 model */
  float boost_beta;

  /*! Normalisation density (internal units) for the boost of the Booth & Schaye
   * 2009 model */
  double boost_n_h_star;

  /*! Switch for nibbling mode */
  int use_nibbling;

  /*! Minimum gas particle mass in nibbling mode */
  float min_gas_mass_for_nibbling;

  /* ---- Properties of the feedback model ------- */

  /*! AGN feedback model: random, isotropic or minimum distance */
  enum AGN_feedback_models feedback_model;

  /*! Is the AGN feedback model deterministic or stochastic? */
  int AGN_deterministic;

  /*! Feedback coupling efficiency of the black holes. */
  float epsilon_f;

  /*! (Constant) temperature increase induced by AGN feedback [Kelvin], if we
   * use a model with a variable temperature increase than we use this value
   * to initialize a BH that just has formed */
  float AGN_delta_T_desired;

  /*! Switch on adaptive heating temperature scheme? */
  int use_variable_delta_T;

  /*! If we use variable delta_T, should we scale with local gas properties
   *  in addition to BH mass? */
  int AGN_with_locally_adaptive_delta_T;

  /*! Normalisation for dT scaling with BH mass */
  float AGN_delta_T_mass_norm;

  /*! Reference BH mass for dT scaling [M_Sun] */
  float AGN_delta_T_mass_reference;

  /*! Exponent for dT scaling with BH mass */
  float AGN_delta_T_mass_exponent;

  /*! Buffer factor for numerical efficiency temperature */
  float AGN_delta_T_crit_factor;

  /*! Buffer factor for background temperature */
  float AGN_delta_T_background_factor;

  /*! Max/min temperature increase induced by AGN feedback [Kelvin] */
  float AGN_delta_T_max;
  float AGN_delta_T_min;

  /*! Vary the energy reservoir according to the BH accretion rate? */
  int use_adaptive_energy_reservoir_threshold;

  /*! Normalisation for energy reservoir threshold, at upper end */
  float nheat_alpha;

  /*! Reference max accretion rate for energy reservoir variation */
  float nheat_maccr_normalisation;

  /*! Hard limit to the energy reservoir threshold */
  float nheat_limit;

  /*! Number of gas neighbours to heat in a feedback event */
  float num_ngbs_to_heat;

  /*! Switch to make nheat use the constant dT as basis, not actual dT */
  int AGN_use_nheat_with_fixed_dT;

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

  /*! Reference black hole mass for repositioning scaling */
  float reposition_reference_mass;

  /*! Repositioning velocity scaling with black hole mass */
  float reposition_exponent_mass;

  /*! Reference gas density for repositioning scaling */
  float reposition_reference_n_H;

  /*! Repositioning velocity scaling with gas density */
  float reposition_exponent_n_H;

  /* ---- Properties of the merger model ---------- */

  /*! Mass ratio above which a merger is considered 'minor' */
  float minor_merger_threshold;

  /*! Mass ratio above which a merger is considered 'major' */
  float major_merger_threshold;

  /*! Type of merger threshold (0: standard, 1: improved) */
  int merger_threshold_type;

  /*! Maximal distance over which BHs merge, in units of softening length */
  float max_merging_distance_ratio;

  /* ---- Black hole time-step properties ---------- */

  /*! Minimum allowed time-step of BH (internal units) */
  float time_step_min;

  /* ---- Common conversion factors --------------- */

  /*! Conversion factor from temperature to internal energy */
  float temp_to_u_factor;

  /*! Conversion factor from physical density to n_H [cgs] */
  float rho_to_n_cgs;

  /*! Conversion factor from internal mass to solar masses */
  float mass_to_solar_mass;
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
      parser_get_param_float(params, "EAGLEAGN:subgrid_seed_mass_Msun");

  /* Convert to internal units */
  bp->subgrid_seed_mass *= phys_const->const_solar_mass;

  bp->use_subgrid_mass_from_ics =
      parser_get_opt_param_int(params, "EAGLEAGN:use_subgrid_mass_from_ics", 1);
  if (bp->use_subgrid_mass_from_ics)
    bp->with_subgrid_mass_check =
        parser_get_opt_param_int(params, "EAGLEAGN:with_subgrid_mass_check", 1);

  /* Accretion parameters ---------------------------------- */

  bp->use_multi_phase_bondi =
      parser_get_param_int(params, "EAGLEAGN:use_multi_phase_bondi");

  bp->use_subgrid_bondi =
      parser_get_param_int(params, "EAGLEAGN:use_subgrid_bondi");

  if (bp->use_multi_phase_bondi && bp->use_subgrid_bondi)
    error(
        "Cannot run with both the multi-phase Bondi and subgrid Bondi models "
        "at the same time!");

  /* Rosas-Guevara et al. (2015) model */
  bp->with_angmom_limiter =
      parser_get_param_int(params, "EAGLEAGN:with_angmom_limiter");
  if (bp->with_angmom_limiter)
    bp->alpha_visc = parser_get_param_float(params, "EAGLEAGN:viscous_alpha");

  bp->epsilon_r =
      parser_get_param_float(params, "EAGLEAGN:radiative_efficiency");
  if (bp->epsilon_r > 1.f)
    error("EAGLEAGN:radiative_efficiency must be <= 1, not %f.", bp->epsilon_r);

  bp->f_Edd = parser_get_param_float(params, "EAGLEAGN:max_eddington_fraction");
  bp->f_Edd_recording = parser_get_param_float(
      params, "EAGLEAGN:eddington_fraction_for_recording");

  /*  Booth & Schaye (2009) Parameters */
  bp->with_boost_factor =
      parser_get_param_int(params, "EAGLEAGN:with_boost_factor");

  if (bp->with_boost_factor) {
    bp->boost_alpha_only =
        parser_get_param_int(params, "EAGLEAGN:boost_alpha_only");
    bp->boost_alpha = parser_get_param_float(params, "EAGLEAGN:boost_alpha");

    if (!bp->boost_alpha_only) {
      bp->boost_beta = parser_get_param_float(params, "EAGLEAGN:boost_beta");

      /* Load the density in cgs and convert to internal units */
      bp->boost_n_h_star =
          parser_get_param_float(params, "EAGLEAGN:boost_n_h_star_H_p_cm3") /
          units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);
    }
  }

  bp->use_nibbling = parser_get_param_int(params, "EAGLEAGN:use_nibbling");
  if (bp->use_nibbling) {
    bp->min_gas_mass_for_nibbling =
        parser_get_param_float(params, "EAGLEAGN:min_gas_mass_for_nibbling");
    bp->min_gas_mass_for_nibbling *= phys_const->const_solar_mass;
  }

  /* Feedback parameters ---------------------------------- */

  char temp[40];
  parser_get_param_string(params, "EAGLEAGN:AGN_feedback_model", temp);
  if (strcmp(temp, "Random") == 0)
    bp->feedback_model = AGN_random_ngb_model;
  else if (strcmp(temp, "Isotropic") == 0)
    bp->feedback_model = AGN_isotropic_model;
  else if (strcmp(temp, "MinimumDistance") == 0)
    bp->feedback_model = AGN_minimum_distance_model;
  else if (strcmp(temp, "MinimumDensity") == 0)
    bp->feedback_model = AGN_minimum_density_model;
  else
    error(
        "The AGN feedback model must be either 'Random', 'MinimumDistance', "
        "'MinimumDensity' or 'Isotropic', not %s",
        temp);

  bp->AGN_deterministic =
      parser_get_param_int(params, "EAGLEAGN:AGN_use_deterministic_feedback");

  bp->epsilon_f =
      parser_get_param_float(params, "EAGLEAGN:coupling_efficiency");

  const double T_K_to_int =
      1. / units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  /* Read the constant AGN heating temperature or the the initial value
   * for the IC or new BH that formed from gas */
  bp->AGN_delta_T_desired =
      parser_get_param_float(params, "EAGLEAGN:AGN_delta_T_K");

  /* Read the properties of the variable heating temperature model */
  bp->use_variable_delta_T =
      parser_get_param_int(params, "EAGLEAGN:use_variable_delta_T");
  if (bp->use_variable_delta_T) {
    bp->AGN_delta_T_mass_norm =
        parser_get_param_float(params, "EAGLEAGN:AGN_delta_T_mass_norm") *
        T_K_to_int;
    bp->AGN_delta_T_mass_reference =
        parser_get_param_float(params, "EAGLEAGN:AGN_delta_T_mass_reference") *
        phys_const->const_solar_mass;
    bp->AGN_delta_T_mass_exponent =
        parser_get_param_float(params, "EAGLEAGN:AGN_delta_T_mass_exponent");

    bp->AGN_with_locally_adaptive_delta_T = parser_get_param_int(
        params, "EAGLEAGN:AGN_with_locally_adaptive_delta_T");
    if (bp->AGN_with_locally_adaptive_delta_T) {
      bp->AGN_delta_T_crit_factor =
          parser_get_param_float(params, "EAGLEAGN:AGN_delta_T_crit_factor");
      bp->AGN_delta_T_background_factor = parser_get_param_float(
          params, "EAGLEAGN:AGN_delta_T_background_factor");
    }

    bp->AGN_delta_T_max =
        parser_get_param_float(params, "EAGLEAGN:AGN_delta_T_max") * T_K_to_int;
    bp->AGN_delta_T_min =
        parser_get_param_float(params, "EAGLEAGN:AGN_delta_T_min") * T_K_to_int;
    bp->AGN_use_nheat_with_fixed_dT =
        parser_get_param_int(params, "EAGLEAGN:AGN_use_nheat_with_fixed_dT");
  }
  bp->use_adaptive_energy_reservoir_threshold = parser_get_param_int(
      params, "EAGLEAGN:AGN_use_adaptive_energy_reservoir_threshold");
  if (bp->use_adaptive_energy_reservoir_threshold) {
    bp->nheat_alpha =
        parser_get_param_float(params, "EAGLEAGN:AGN_nheat_alpha");
    bp->nheat_maccr_normalisation =
        parser_get_param_float(params,
                               "EAGLEAGN:AGN_nheat_maccr_normalisation") *
        phys_const->const_solar_mass / phys_const->const_year;
    bp->nheat_limit =
        parser_get_param_float(params, "EAGLEAGN:AGN_nheat_limit");
  }

  /* We must always read a default value to initialize BHs to */
  bp->num_ngbs_to_heat =
      parser_get_param_float(params, "EAGLEAGN:AGN_num_ngb_to_heat");

  /* Reposition parameters --------------------------------- */

  bp->max_reposition_mass =
      parser_get_param_float(params, "EAGLEAGN:max_reposition_mass") *
      phys_const->const_solar_mass;
  bp->max_reposition_distance_ratio =
      parser_get_param_float(params, "EAGLEAGN:max_reposition_distance_ratio");

  bp->with_reposition_velocity_threshold = parser_get_param_int(
      params, "EAGLEAGN:with_reposition_velocity_threshold");

  if (bp->with_reposition_velocity_threshold) {
    bp->max_reposition_velocity_ratio = parser_get_param_float(
        params, "EAGLEAGN:max_reposition_velocity_ratio");

    /* Prevent nonsensical input */
    if (bp->max_reposition_velocity_ratio <= 0)
      error("max_reposition_velocity_ratio must be positive, not %f.",
            bp->max_reposition_velocity_ratio);

    bp->min_reposition_velocity_threshold = parser_get_param_float(
        params, "EAGLEAGN:min_reposition_velocity_threshold");
    /* Convert from km/s to internal units */
    bp->min_reposition_velocity_threshold *=
        (1e5 / (us->UnitLength_in_cgs / us->UnitTime_in_cgs));
  }

  bp->set_reposition_speed =
      parser_get_param_int(params, "EAGLEAGN:set_reposition_speed");

  if (bp->set_reposition_speed) {
    bp->reposition_coefficient_upsilon = parser_get_param_float(
        params, "EAGLEAGN:reposition_coefficient_upsilon");

    /* Prevent the user from making silly wishes */
    if (bp->reposition_coefficient_upsilon <= 0)
      error(
          "reposition_coefficient_upsilon must be positive, not %f "
          "km/s/M_sun.",
          bp->reposition_coefficient_upsilon);

    /* Convert from km/s to internal units */
    bp->reposition_coefficient_upsilon *=
        (1e5 / (us->UnitLength_in_cgs / us->UnitTime_in_cgs));

    /* Scaling parameters with BH mass and gas density */
    bp->reposition_reference_mass =
        parser_get_param_float(params, "EAGLEAGN:reposition_reference_mass") *
        phys_const->const_solar_mass;
    bp->reposition_exponent_mass = parser_get_opt_param_float(
        params, "EAGLEAGN:reposition_exponent_mass", 2.0);
    bp->reposition_reference_n_H =
        parser_get_param_float(params, "EAGLEAGN:reposition_reference_n_H");
    bp->reposition_exponent_n_H = parser_get_opt_param_float(
        params, "EAGLEAGN:reposition_exponent_n_H", 1.0);
  }

  /* Merger parameters ------------------------------------- */

  bp->minor_merger_threshold =
      parser_get_param_float(params, "EAGLEAGN:threshold_minor_merger");

  bp->major_merger_threshold =
      parser_get_param_float(params, "EAGLEAGN:threshold_major_merger");

  bp->merger_threshold_type =
      parser_get_param_int(params, "EAGLEAGN:merger_threshold_type");

  bp->max_merging_distance_ratio =
      parser_get_param_float(params, "EAGLEAGN:merger_max_distance_ratio");

  /* ---- Black hole time-step properties ------------------ */

  const double Myr_in_cgs = 1e6 * 365.25 * 24. * 60. * 60.;

  const double time_step_min_Myr = parser_get_opt_param_float(
      params, "EAGLEAGN:minimum_timestep_Myr", FLT_MAX);

  bp->time_step_min = time_step_min_Myr * Myr_in_cgs /
                      units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /* Common conversion factors ----------------------------- */

  /* Calculate temperature to internal energy conversion factor (all internal
   * units) */
  const double k_B = phys_const->const_boltzmann_k;
  const double m_p = phys_const->const_proton_mass;
  const double mu = hydro_props->mu_ionised;
  bp->temp_to_u_factor = k_B / (mu * hydro_gamma_minus_one * m_p);

  /* Calculate conversion factor from rho to n_H.
   * Note this assumes primoridal abundance */
  const double X_H = hydro_props->hydrogen_mass_fraction;
  bp->rho_to_n_cgs =
      (X_H / m_p) * units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);

  /* Conversion factor for internal mass to M_solar */
  bp->mass_to_solar_mass = 1. / phys_const->const_solar_mass;
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

#endif /* SWIFT_EAGLE_BLACK_HOLES_PROPERTIES_H */
