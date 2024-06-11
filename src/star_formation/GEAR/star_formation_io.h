/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
 *               2019 Fabien Jeanquartier (fabien.jeanquartier@epfl.ch)
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
#ifndef SWIFT_STAR_FORMATION_GEAR_IO_H
#define SWIFT_STAR_FORMATION_GEAR_IO_H

/* Config parameters. */
#include <config.h>

/* Local includes */
#include "io_properties.h"

/**
 * @brief Specifies which s-particle fields to read from a dataset
 *
 * @param sparts The s-particle array.
 * @param list The list of i/o properties to read.
 *
 * @return num_fields The number of i/o fields to read.
 */
INLINE static int star_formation_read_particles(struct spart* sparts,
                                                struct io_props* list) {

  /* List what we want to read */
  list[0] = io_make_input_field("BirthMass", FLOAT, 1, OPTIONAL, UNIT_CONV_MASS,
                                sparts, sf_data.birth_mass);

  return 1;
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extended data particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
__attribute__((always_inline)) INLINE static int star_formation_write_particles(
    const struct part* parts, const struct xpart* xparts,
    struct io_props* list) {
  /* Nothing to write here */
  return 0;
}

/**
 * @brief Specifies which sparticle fields to write to a dataset
 *
 * @param sparts The star particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
__attribute__((always_inline)) INLINE static int
star_formation_write_sparticles(const struct spart* sparts,
                                struct io_props* list) {

  list[0] = io_make_output_field(
      "BirthDensities", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, sparts,
      sf_data.birth_density,
      "Physical densities at the time of birth of the gas particles that "
      "turned into stars (note that "
      "we store the physical density at the birth redshift, no conversion is "
      "needed)");

  list[1] =
      io_make_output_field("BirthTemperatures", FLOAT, 1, UNIT_CONV_TEMPERATURE,
                           0.f, sparts, sf_data.birth_temperature,
                           "Temperatures at the time of birth of the gas "
                           "particles that turned into stars");

  list[2] = io_make_output_field("BirthMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f,
                                 sparts, sf_data.birth_mass,
                                 "Masses of the star particles at birth time");

  list[3] = io_make_output_field(
      "ProgenitorIDs", LONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      sf_data.progenitor_id, "Unique IDs of the progenitor particle");

  list[4] =
      io_make_output_field("StellarParticleType", CHAR, 1, UNIT_CONV_NO_UNITS,
                           0.f, sparts, feedback_data.star_type,
                           "Type of stellar particle: 0=single star ; 1=stellar"
                           " part. without SNII ; 2=normal");

  return 5;
}

/**
 * @brief initialization of the star formation law
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param hydro_props The #hydro_props.
 * @param cosmo The current cosmological model.
 * @param entropy_floor The properties of the entropy floor used in this
 * simulation.
 * @param starform the star formation law properties to initialize
 */
INLINE static void starformation_init_backend(
    struct swift_params* parameter_file, const struct phys_const* phys_const,
    const struct unit_system* us, const struct hydro_props* hydro_props,
    const struct cosmology* cosmo,
    const struct entropy_floor_properties* entropy_floor,
    struct star_formation* starform) {

  const char* default_mode = "default";
  char temp[32];
  parser_get_opt_param_string(parameter_file,
                              "GEARStarFormation:star_formation_mode", temp,
                              default_mode);

  /* Star formation mode */
  if (strcmp(temp, "default") == 0) {
    starform->star_formation_mode = gear_star_formation_default;
  } else if (strcmp(temp, "agora") == 0) {
    starform->star_formation_mode = gear_star_formation_agora;
  } else {
    error("Invalid star formation model: '%s'", temp);
  }

  /* Star formation efficiency */
  starform->star_formation_efficiency = parser_get_param_double(
      parameter_file, "GEARStarFormation:star_formation_efficiency");

  /* Maximum gas temperature for star formation */
  starform->maximal_temperature = parser_get_param_double(
      parameter_file, "GEARStarFormation:maximal_temperature");

  /* Minimal gas density for star formation */
  starform->density_threshold = parser_get_param_double(
      parameter_file, "GEARStarFormation:density_threshold");

  /* Number of stars per particles */
  starform->n_stars_per_part = parser_get_param_double(
      parameter_file, "GEARStarFormation:n_stars_per_particle");

  /* Minimal fraction of mass for the last star formed. */
  starform->min_mass_frac_plus_one = parser_get_param_double(
      parameter_file, "GEARStarFormation:min_mass_frac");
  /* Avoid generating gas particle with mass below the fraction => + 1. */
  starform->min_mass_frac_plus_one += 1.;

  /* Get the jeans factor */
  starform->n_jeans_2_3 =
      parser_get_param_float(parameter_file, "GEARPressureFloor:jeans_factor");
  starform->n_jeans_2_3 = pow(starform->n_jeans_2_3, 2. / 3.);

  /* Apply unit change */
  starform->maximal_temperature /=
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  starform->density_threshold /=
      units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);

  /* Initialize the mass of the stars to 0 for the stats computation */
  starform->mass_stars = 0;

  /* Print parameters */
  if (engine_rank == 0) {
    message("star_formation_mode       = %d", starform->star_formation_mode);
    message("star_formation_efficiency = %g",
            starform->star_formation_efficiency);
    message("maximal_temperature       = %g", starform->maximal_temperature);
    message("density_threshold         = %g", starform->density_threshold);
    message("n_stars_per_part          = %d", starform->n_stars_per_part);
    message("min_mass_frac_plus_one    = %g", starform->min_mass_frac_plus_one);
    message("n_jeans_2_3               = %g", starform->n_jeans_2_3);
  }
}

#endif /* SWIFT_STAR_FORMATION_GEAR_IO_H */
