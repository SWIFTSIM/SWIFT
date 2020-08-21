/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_STAR_FORMATION_EAGLE_IO_H
#define SWIFT_STAR_FORMATION_EAGLE_IO_H

/* Config parameters. */
#include "../config.h"

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
  return 0;
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

  list[0] = io_make_output_field(
      "StarFormationRates", FLOAT, 1, UNIT_CONV_SFR, 0.f, xparts, sf_data.SFR,
      "If positive, star formation rates of the particles. If negative, stores "
      "the last time/scale-factor at which the gas particle was star-forming. "
      "If zero, the particle was never star-forming.");

  return 1;
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

  return 2;
}
#endif /* SWIFT_STAR_FORMATION_EAGLE_IO_H */
