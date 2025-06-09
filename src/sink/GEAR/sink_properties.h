/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
 *               2024 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_GEAR_SINK_PROPERTIES_H
#define SWIFT_GEAR_SINK_PROPERTIES_H

/* Local header */
#include "feedback_properties.h"
#include "parser.h"

/* Some default values for the parameters to be read in the YAML file */
#define sink_gear_f_acc_default 0.8
#define sink_gear_star_spawning_sigma_factor_default 0.2
#define sink_gear_n_imf_default FLT_MAX /* No accretion restriction */
#define sink_gear_tolerance_sf_timestep_default 0.5

/* Sink formation is activated */
#define sink_gear_disable_sink_formation_default 0

/* By default all current implemented criteria are active */
#define sink_gear_sink_formation_criterion_all_default 1

/**
 * @brief Properties of sink in the Default model.
 */
struct sink_props {

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

  /*! Are we using a fixed cutoff radius? (all smoothing length calculations are
   * disabled if so) */
  char use_fixed_r_cut;

  /*! Cut off radius */
  float cut_off_radius;

  /* Fraction of the cut-off radius for gas accretion. It should respect 0 <=
     f_acc <= 1 */
  float f_acc;

  /*! Maximal gas temperature for forming a sink. */
  float temperature_threshold;

  /*! Minimal gas density for forming a sink with the temperature threshold. */
  float density_threshold;

  /*! Gas density for forming a sink without the temperature threshold. */
  float maximal_density_threshold;

  /*! Mass of the stellar particle representing the low mass stars
   * (continuous IMF sampling). In M_sun. */
  float stellar_particle_mass_Msun;

  /*! Minimal mass of stars represented by discrete particles. In M_sun. */
  float minimal_discrete_mass_Msun;

  /*! Mass of the stellar particle representing the low mass stars
   * (continuous IMF sampling). In M_sun. First stars */
  float stellar_particle_mass_first_stars_Msun;

  /*! Minimal mass of stars represented by discrete particles. In M_sun.
   * First stars. */
  float minimal_discrete_mass_first_stars_Msun;

  /*! Sink formation criteria selecter : some criteria can be left out.  */
  char sink_formation_contracting_gas_criterion;
  char sink_formation_smoothing_length_criterion;
  char sink_formation_jeans_instability_criterion;
  char sink_formation_bound_state_criterion;
  char sink_formation_overlapping_sink_criterion;

  /*! Disable sink formation? (e.g. used in sink accretion tests). Default: 0
     (keep sink formation) */
  char disable_sink_formation;

  /*! Factor to rescale the velocity dispersion of the stars when they are
     spawned */
  double star_spawning_sigma_factor;

  /*! Minimal sink mass in Msun. This prevents m_sink << m_gas in low
    resolution simulations. */
  float sink_minimal_mass_Msun;

  /***************************************************************************/
  /*! Maximal time-step length of young sinks (internal units) */
  double max_time_step_young;

  /*! Maximal time-step length of old sinks (internal units) */
  double max_time_step_old;

  /*! Age threshold for the young/old transition (internal units) */
  double age_threshold;

  /*! Age threshold for the transition to unlimited time-step size (internal
   * units) */
  double age_threshold_unlimited;

  /*! Time integration CFL condition factor */
  float CFL_condition;

  /*! Number of times the IMF mass can be swallowed in a single timestep */
  float n_IMF;

  /*! Tolerance parameter for SF timestep constraint */
  float tolerance_SF_timestep;
};

/**
 * @brief Initialise the probabilities to get a stellar mass (continuous
 * sampling of the IMF)
 *
 * @param sp The #sink_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param cosmo The cosmological model.
 */
INLINE static void sink_props_init_probabilities(
    struct sink_props *sp, struct initial_mass_function *imf,
    const struct phys_const *phys_const, int first_stars) {

  /* get the IMF mass limits (all in Msol) */
  float mass_min = imf->mass_min;
  float mass_max = imf->mass_max;

  /* Treat separately the cases of first star or not. */
  if (!first_stars) {
    imf->minimal_discrete_mass_Msun = sp->minimal_discrete_mass_Msun;
    imf->stellar_particle_mass_Msun = sp->stellar_particle_mass_Msun;
  } else {
    imf->minimal_discrete_mass_Msun =
        sp->minimal_discrete_mass_first_stars_Msun;
    imf->stellar_particle_mass_Msun =
        sp->stellar_particle_mass_first_stars_Msun;
  }

  /* Sanity check */
  if (imf->minimal_discrete_mass_Msun < imf->mass_limits[imf->n_parts - 1])
    error(
        "minimal_discrete_mass (=%8.3f) cannot be smaller than the mass limit "
        "(=%8.3f) of the last IMF segment,",
        imf->minimal_discrete_mass_Msun, imf->mass_limits[imf->n_parts - 1]);

  /* Compute the IMF mass (in solar mass) below the minimal IMF discrete mass
     (continuous part). */
  double Mtot, Md, Mc;
  initial_mass_function_compute_Mc_Md_Mtot(imf, &Mc, &Md, &Mtot);

  /* Compute the number of stars in the continuous part of the IMF */
  double Nc = initial_mass_function_get_imf_number_fraction(
                  imf, mass_min, imf->minimal_discrete_mass_Msun) *
              Mtot;

  /* Compute the number of stars in the discrete part of the IMF */
  double Nd = initial_mass_function_get_imf_number_fraction(
                  imf, imf->minimal_discrete_mass_Msun, mass_max) *
              Mtot;

  if (engine_rank == 0) {
    message("Mass of the continuous part (in M_sun) : %g", Mc);
    message("Mass of the discrete   part (in M_sun) : %g", Md);
    message("Total IMF mass (in M_sun)              : %g", Mtot);
    message("Number of stars in the continuous part : %g", Nc);
    message("Number of stars in the discrete   part : %g", Nd);
  }

  /* if no continous part, return */
  if (Mc == 0) {
    imf->sink_Pc = 0;
    imf->stellar_particle_mass_Msun = 0;
    if (engine_rank == 0) {
      message("probability of the continuous part    : %g", 0.);
      message("probability of the discrete   part    : %g", 1.);
    }
    return;
  }

  /* Compute the probabilities */
  double Pc = 1 / (1 + Nd);
  double Pd = 1 - Pc;
  imf->sink_Pc = Pc;

  if (engine_rank == 0) {
    message("probability of the continuous part     : %g", Pc);
    message("probability of the discrete   part     : %g", Pd);
  }
}

/**
 * @brief Initialise the sink properties from the parameter file.
 *
 * @param sp The #sink_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param cosmo The cosmological model.
 * @param with_feedback Are we running with feedback?
 */
INLINE static void sink_props_init(
    struct sink_props *sp, struct feedback_props *fp,
    const struct phys_const *phys_const, const struct unit_system *us,
    struct swift_params *params, const struct hydro_props *hydro_props,
    const struct cosmology *cosmo, const int with_feedback) {

  /* Read in the basic neighbour search properties or default to the hydro
     ones if the user did not provide any different values */

  /* Kernel properties */
  sp->eta_neighbours = parser_get_opt_param_float(
      params, "Sinks:resolution_eta", hydro_props->eta_neighbours);

  /* Tolerance for the smoothing length Newton-Raphson scheme */
  sp->h_tolerance = parser_get_opt_param_float(params, "Sinks:h_tolerance",
                                               hydro_props->h_tolerance);

  /* Get derived properties */
  sp->target_neighbours = pow_dimension(sp->eta_neighbours) * kernel_norm;
  const float delta_eta = sp->eta_neighbours * (1.f + sp->h_tolerance);
  sp->delta_neighbours =
      (pow_dimension(delta_eta) - pow_dimension(sp->eta_neighbours)) *
      kernel_norm;

  /* Number of iterations to converge h */
  sp->max_smoothing_iterations =
      parser_get_opt_param_int(params, "Sinks:max_ghost_iterations",
                               hydro_props->max_smoothing_iterations);

  /* Time integration properties */
  const float max_volume_change =
      parser_get_opt_param_float(params, "Sinks:max_volume_change", -1);
  if (max_volume_change == -1)
    sp->log_max_h_change = hydro_props->log_max_h_change;
  else
    sp->log_max_h_change = logf(powf(max_volume_change, hydro_dimension_inv));

  /* If we do not run with feedback, abort and print an error */
  if (!with_feedback)
    error(
        "ERROR: Running with sink but without feedback. GEAR sink model needs "
        "to be run with --sink and --feedback");

  /* This property is used in all models to flag if we're using a fixed cutoff.
   * If it is set to 1, we use a fixed r_cut (read in below) and don't
   * (re)calculate h. If not, the code will use the variable h*kernel_gamma as a
   * cutoff radius.
   */
  sp->use_fixed_r_cut =
      parser_get_param_char(params, "GEARSink:use_fixed_cut_off_radius");

  /* The property cut_off_radius is now only used in the GEAR model.
   * It is ignored if use_fixed_r_cut is 0. */
  if (sp->use_fixed_r_cut) {
    sp->cut_off_radius =
        parser_get_param_float(params, "GEARSink:cut_off_radius");
  }

  sp->f_acc = parser_get_opt_param_float(params, "GEARSink:f_acc",
                                         sink_gear_f_acc_default);

  /* Check that sp->f_acc respects 0 <= f_acc <= 1 */
  if ((sp->f_acc < 0) || (sp->f_acc > 1)) {
    error(
        "The sink f_acc has not an allowed value. It should respect 0 <= f_acc "
        "<= 1. Current value f_acc = %f.",
        sp->f_acc);
  }

  sp->temperature_threshold =
      parser_get_param_float(params, "GEARSink:temperature_threshold_K");

  sp->density_threshold =
      parser_get_param_float(params, "GEARSink:density_threshold_Hpcm3");

  sp->maximal_density_threshold = parser_get_opt_param_float(
      params, "GEARSink:maximal_density_threshold_Hpcm3", FLT_MAX);

  if (sp->maximal_density_threshold < sp->density_threshold) {
    error(
        "maximal_density_threshold_Hpcm3 must be larger than "
        "density_threshold_Hpcm3");
  }

  sp->stellar_particle_mass_Msun =
      parser_get_param_float(params, "GEARSink:stellar_particle_mass_Msun");

  sp->minimal_discrete_mass_Msun =
      parser_get_param_float(params, "GEARSink:minimal_discrete_mass_Msun");

  sp->stellar_particle_mass_first_stars_Msun = parser_get_param_float(
      params, "GEARSink:stellar_particle_mass_first_stars_Msun");

  sp->minimal_discrete_mass_first_stars_Msun = parser_get_param_float(
      params, "GEARSink:minimal_discrete_mass_first_stars_Msun");

  sp->star_spawning_sigma_factor =
      parser_get_opt_param_float(params, "GEARSink:star_spawning_sigma_factor",
                                 sink_gear_star_spawning_sigma_factor_default);

  sp->n_IMF = parser_get_opt_param_float(params, "GEARSink:n_IMF",
                                         sink_gear_n_imf_default);

  sp->sink_minimal_mass_Msun =
      parser_get_opt_param_float(params, "GEARSink:sink_minimal_mass_Msun", 0.);

  /* Sink formation criterion parameters (all active by default) */
  sp->sink_formation_contracting_gas_criterion = parser_get_opt_param_int(
      params, "GEARSink:sink_formation_contracting_gas_criterion",
      sink_gear_sink_formation_criterion_all_default);

  sp->sink_formation_smoothing_length_criterion = parser_get_opt_param_int(
      params, "GEARSink:sink_formation_smoothing_length_criterion",
      sink_gear_sink_formation_criterion_all_default);

  sp->sink_formation_jeans_instability_criterion = parser_get_opt_param_int(
      params, "GEARSink:sink_formation_jeans_instability_criterion",
      sink_gear_sink_formation_criterion_all_default);

  sp->sink_formation_bound_state_criterion = parser_get_opt_param_int(
      params, "GEARSink:sink_formation_bound_state_criterion",
      sink_gear_sink_formation_criterion_all_default);

  sp->sink_formation_overlapping_sink_criterion = parser_get_opt_param_int(
      params, "GEARSink:sink_formation_overlapping_sink_criterion",
      sink_gear_sink_formation_criterion_all_default);

  /* Should we disable sink formation ? */
  sp->disable_sink_formation =
      parser_get_opt_param_int(params, "GEARSink:disable_sink_formation",
                               sink_gear_disable_sink_formation_default);

  /* Maximal time-step lengths */
  const double max_time_step_young_Myr = parser_get_opt_param_float(
      params, "GEARSink:max_timestep_young_Myr", FLT_MAX);
  const double max_time_step_old_Myr = parser_get_opt_param_float(
      params, "GEARSink:max_timestep_old_Myr", FLT_MAX);
  const double age_threshold_Myr = parser_get_opt_param_float(
      params, "GEARSink:timestep_age_threshold_Myr", FLT_MAX);
  const double age_threshold_unlimited_Myr = parser_get_opt_param_float(
      params, "GEARSink:timestep_age_threshold_unlimited_Myr", FLT_MAX);

  /* Check for consistency */
  if (age_threshold_unlimited_Myr != 0. && age_threshold_Myr != FLT_MAX) {
    if (age_threshold_unlimited_Myr < age_threshold_Myr)
      error(
          "The age threshold for unlimited sink time-step sizes (%e Myr) is "
          "smaller than the transition threshold from young to old ages (%e "
          "Myr)",
          age_threshold_unlimited_Myr, age_threshold_Myr);
  }

  /* Timestep tolerance paramters */
  sp->CFL_condition = parser_get_param_float(params, "GEARSink:CFL_condition");

  sp->tolerance_SF_timestep =
      parser_get_opt_param_float(params, "GEARSink:tolerance_SF_timestep",
                                 sink_gear_tolerance_sf_timestep_default);

  /* Apply unit change */
  sp->temperature_threshold /=
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  const double m_p_cgs = phys_const->const_proton_mass *
                         units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  sp->density_threshold *=
      m_p_cgs / units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  sp->maximal_density_threshold *=
      m_p_cgs / units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);

  const double Myr_internal_units = 1e6 * phys_const->const_year;
  sp->max_time_step_young = max_time_step_young_Myr * Myr_internal_units;
  sp->max_time_step_old = max_time_step_old_Myr * Myr_internal_units;
  sp->age_threshold = age_threshold_Myr * Myr_internal_units;
  sp->age_threshold_unlimited =
      age_threshold_unlimited_Myr * Myr_internal_units;

  /* here, we need to differenciate between the stellar models */
  struct initial_mass_function *imf;
  struct stellar_model *sm;

  sm = &fp->stellar_model;
  imf = &sm->imf;

  /* Initialize for the stellar models (PopII) */
  sink_props_init_probabilities(sp, imf, phys_const, 0);

  /* Now initialize the first stars. */
  if (fp->metallicity_max_first_stars != -1) {
    sm = &fp->stellar_model_first_stars;
    imf = &sm->imf;
    sink_props_init_probabilities(sp, imf, phys_const, 1);
  }
  if (engine_rank == 0) {
    message("temperature_threshold                        = %g",
            sp->temperature_threshold);
    message("density_threshold                            = %g",
            sp->density_threshold);
    message("maximal_density_threshold                    = %g",
            sp->maximal_density_threshold);
    message("sink_minimal_mass (in M_sun)                 = %g",
            sp->sink_minimal_mass_Msun);
    message("stellar_particle_mass (in M_sun)             = %g",
            sp->stellar_particle_mass_Msun);
    message("minimal_discrete_mass (in M_sun)             = %g",
            sp->minimal_discrete_mass_Msun);

    message("stellar_particle_mass_first_stars (in M_sun) = %g",
            sp->stellar_particle_mass_first_stars_Msun);
    message("minimal_discrete_mass_first_stars (in M_sun) = %g",
            sp->minimal_discrete_mass_first_stars_Msun);

    /* Print information about the functionalities */
    message("disable_sink_formation                       = %d",
            sp->disable_sink_formation);
    message("sink_formation_contracting_gas_criterion     = %d",
            sp->sink_formation_contracting_gas_criterion);
    message("sink_formation_smoothing_length_criterion    = %d",
            sp->sink_formation_smoothing_length_criterion);
    message("sink_formation_jeans_instability_criterion   = %d",
            sp->sink_formation_jeans_instability_criterion);
    message("sink_formation_bound_state_criterion         = %d",
            sp->sink_formation_bound_state_criterion);
    message("sink_formation_overlapping_sink_criterion    = %d",
            sp->sink_formation_overlapping_sink_criterion);

    /* Print timestep parameters information */
    message("sink max_timestep_young                      = %e",
            sp->max_time_step_young);
    message("sink max_timestep_old                        = %e",
            sp->max_time_step_old);
    message("sink age_threshold from young to old         = %e",
            sp->age_threshold);
    message("sink age_threshold from old to unlimited     = %e",
            sp->age_threshold_unlimited);
    message("sink C_CFL                                   = %e",
            sp->CFL_condition);
    message("tolerance_SF_timestep                        = %e",
            sp->tolerance_SF_timestep);
    message("n_IMF                                        = %e", sp->n_IMF);
  }
}

/**
 * @brief Write a sink_props struct to the given FILE as a stream of
 * bytes.
 *
 * @param props the sink properties struct
 * @param stream the file stream
 */
INLINE static void sink_struct_dump(const struct sink_props *props,
                                    FILE *stream) {
  restart_write_blocks((void *)props, sizeof(struct sink_props), 1, stream,
                       "sink props", "Sink props");
}

/**
 * @brief Restore a sink_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param props the sink properties struct
 * @param stream the file stream
 */
INLINE static void sink_struct_restore(const struct sink_props *props,
                                       FILE *stream) {
  restart_read_blocks((void *)props, sizeof(struct sink_props), 1, stream, NULL,
                      "Sink props");
}

#endif /* SWIFT_GEAR_SINK_PROPERTIES_H */
