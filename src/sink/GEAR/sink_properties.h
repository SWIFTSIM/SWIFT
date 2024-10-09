/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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

/**
 * @brief Properties of sink in the Default model.
 */
struct sink_props {

  /*! Cut off radius */
  float cut_off_radius;

  /* Fraction of the cut-off radius for gas accretion. It should respect 0 <=
     f_acc <= 1 */
  float f_acc;

  /*! Maximal gas temperature for forming a star. */
  float maximal_temperature;

  /*! Minimal gas density for forming a star. */
  float density_threshold;

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

  /* Disable sink formation? (e.g. used in sink accretion tests). Default: 0
     (keep sink formation) */
  uint8_t disable_sink_formation;

  /* Factor to rescale the velocity dispersion of the stars when they are
     spawned */
  double star_spawning_sigma_factor;
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

  message("Mass of the continuous part (in M_sun) : %g", Mc);
  message("Mass of the discrete   part (in M_sun) : %g", Md);
  message("Total IMF mass (in M_sun)              : %g", Mtot);
  message("Number of stars in the continuous part : %g", Nc);
  message("Number of stars in the discrete   part : %g", Nd);

  /* if no continous part, return */
  if (Mc == 0) {
    imf->sink_Pc = 0;
    imf->stellar_particle_mass_Msun = 0;
    message("probability of the continuous part    : %g", 0.);
    message("probability of the discrete   part    : %g", 1.);
    return;
  }

  /* Compute the probabilities */
  double Pc = 1 / (1 + Nd);
  double Pd = 1 - Pc;
  imf->sink_Pc = Pc;

  message("probability of the continuous part     : %g", Pc);
  message("probability of the discrete   part     : %g", Pd);
}

/**
 * @brief Initialise the sink properties from the parameter file.
 *
 * @param sp The #sink_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param cosmo The cosmological model.
 */
INLINE static void sink_props_init(struct sink_props *sp,
                                   struct feedback_props *fp,
                                   const struct phys_const *phys_const,
                                   const struct unit_system *us,
                                   struct swift_params *params,
                                   const struct cosmology *cosmo) {

  /* Default values */
  const float default_f_acc = 0.8;
  const float default_star_spawning_sigma_factor = 0.2;
  const char default_disable_sink_formation = 0; /* Sink formation is
                                                     activated */

  /* By default all current implemented criteria are active */
  const uint8_t default_sink_formation_criterion_all = 1;

  /* Read the parameters from the parameter file */
  sp->cut_off_radius =
      parser_get_param_float(params, "GEARSink:cut_off_radius");

  sp->f_acc =
      parser_get_opt_param_float(params, "GEARSink:f_acc", default_f_acc);

  /* Check that sp->f_acc respects 0 <= f_acc <= 1 */
  if ((sp->f_acc < 0) || (sp->f_acc > 1)) {
    error(
        "The sink f_acc has not an allowed value. It should respect 0 <= f_acc "
        "<= 1. Current value f_acc = %f.",
        sp->f_acc);
  }

  sp->maximal_temperature =
      parser_get_param_float(params, "GEARSink:maximal_temperature");

  sp->density_threshold =
      parser_get_param_float(params, "GEARSink:density_threshold_g_per_cm3");

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
                                 default_star_spawning_sigma_factor);

  /* Sink formation criterion parameters (all active by default) */
  sp->sink_formation_contracting_gas_criterion = parser_get_opt_param_int(
      params, "GEARSink:sink_formation_contracting_gas_criterion",
      default_sink_formation_criterion_all);

  sp->sink_formation_smoothing_length_criterion = parser_get_opt_param_int(
      params, "GEARSink:sink_formation_smoothing_length_criterion",
      default_sink_formation_criterion_all);

  sp->sink_formation_jeans_instability_criterion = parser_get_opt_param_int(
      params, "GEARSink:sink_formation_jeans_instability_criterion",
      default_sink_formation_criterion_all);

  sp->sink_formation_bound_state_criterion = parser_get_opt_param_int(
      params, "GEARSink:sink_formation_bound_state_criterion",
      default_sink_formation_criterion_all);

  sp->sink_formation_overlapping_sink_criterion = parser_get_opt_param_int(
      params, "GEARSink:sink_formation_overlapping_sink_criterion",
      default_sink_formation_criterion_all);

  /* Should we disable sink formation ? */
  sp->disable_sink_formation =
      parser_get_opt_param_int(params, "GEARSink:disable_sink_formation",
                               default_disable_sink_formation);

  /* Apply unit change */
  sp->maximal_temperature /=
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  sp->density_threshold /= units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);

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

  message("maximal_temperature                          = %g",
          sp->maximal_temperature);
  message("density_threshold                            = %g",
          sp->density_threshold);

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
