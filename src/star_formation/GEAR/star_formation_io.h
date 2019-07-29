/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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
#include "../config.h"

/* Local includes */
#include "io_properties.h"

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
 * @brief initialization of the star formation law
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param starform the star formation law properties to initialize
 *
 */
INLINE static void starformation_init_backend(
    struct swift_params* parameter_file, const struct phys_const* phys_const,
    const struct unit_system* us, const struct hydro_props* hydro_props,
    struct star_formation* starform) {

  // TODO move into pressure floor
  starform->n_jeans_2_3 =
      parser_get_param_float(parameter_file, "GEARStarFormation:NJeans");
  starform->n_jeans_2_3 = pow(starform->n_jeans_2_3, 2. / 3.);

  /* Star formation efficiency */
  starform->star_formation_efficiency = parser_get_param_double(
      parameter_file, "GEARStarFormation:star_formation_efficiency");

  /* Maximum temperature for star formation */
  starform->maximal_temperature = parser_get_param_double(
      parameter_file, "GEARStarFormation:maximal_temperature");

  /* Apply unit change */
  starform->maximal_temperature *=
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
}

#endif /* SWIFT_STAR_FORMATION_GEAR_IO_H */
