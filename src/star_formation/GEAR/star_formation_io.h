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
	/*! we write the properties of our gas particles, star formation rate, temperature, star formation probability, turbulence estimation and the three components of velocity.*/
  list[0] =
      io_make_output_field("SFR", DOUBLE, 1, UNIT_CONV_SFR, xparts, sf_data.SFR);
      list[1]=io_make_output_field("gas_temperature", DOUBLE, 1, UNIT_CONV_TEMPERATURE, xparts, sf_data.temperature);
	  list[2] = io_make_output_field("prob", DOUBLE, 1, UNIT_CONV_NO_UNITS, xparts, sf_data.proba);
      list[3]=io_make_output_field("sigma", FLOAT, 1, UNIT_CONV_VELOCITY, parts, starform_data.sigma_save);
      list[4]=io_make_output_field("velocity1", FLOAT, 1, UNIT_CONV_VELOCITY, parts, v[0]);
      list[5]=io_make_output_field("velocity2", FLOAT, 1, UNIT_CONV_VELOCITY, parts, v[1]);
      list[6]=io_make_output_field("velocity3", FLOAT, 1, UNIT_CONV_VELOCITY, parts, v[2]);
      list[7]=io_make_output_field("density2", DOUBLE, 1, UNIT_CONV_DENSITY, xparts, sf_data.density);
  return 8;
}

#endif /* SWIFT_STAR_FORMATION_GEAR_IO_H */
