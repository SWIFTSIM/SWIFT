/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_TRACERS_EAGLE_IO_H
#define SWIFT_TRACERS_EAGLE_IO_H

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "io_properties.h"
#include "tracers.h"

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of tracers to the file.
 *
 * @param h_grp The HDF5 group in which to write
 */
__attribute__((always_inline)) INLINE static void tracers_write_flavour(
    hid_t h_grp) {

  io_write_attribute_s(h_grp, "Tracers", "EAGLE");
}
#endif

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extended data particle array.
 * @param list The list of i/o properties to write.
 * @param with_cosmology Are we running with cosmology switched on?
 *
 * @return Returns the number of fields to write.
 */
__attribute__((always_inline)) INLINE static int tracers_write_particles(
    const struct part* parts, const struct xpart* xparts, struct io_props* list,
    const int with_cosmology) {

  list[0] = io_make_output_field("Maximal Temperature", FLOAT, 1,
                                 UNIT_CONV_TEMPERATURE, xparts,
                                 tracers_data.maximum_temperature);

  if (with_cosmology) {
    list[1] = io_make_output_field(
        "Maximal Temperature scale-factor", FLOAT, 1, UNIT_CONV_NO_UNITS,
        xparts, tracers_data.maximum_temperature_scale_factor);

  } else {

    list[1] = io_make_output_field("Maximal Temperature time", FLOAT, 1,
                                   UNIT_CONV_NO_UNITS, xparts,
                                   tracers_data.maximum_temperature_time);
  }

  return 2;
}

__attribute__((always_inline)) INLINE static int tracers_write_sparticles(
    const struct spart* sparts, struct io_props* list,
    const int with_cosmology) {

  list[0] = io_make_output_field("Maximal Temperature", FLOAT, 1,
                                 UNIT_CONV_TEMPERATURE, sparts,
                                 tracers_data.maximum_temperature);

  if (with_cosmology) {
    list[1] = io_make_output_field(
        "Maximal Temperature scale-factor", FLOAT, 1, UNIT_CONV_NO_UNITS,
        sparts, tracers_data.maximum_temperature_scale_factor);

  } else {

    list[1] = io_make_output_field("Maximal Temperature time", FLOAT, 1,
                                   UNIT_CONV_NO_UNITS, sparts,
                                   tracers_data.maximum_temperature_time);
  }

  return 2;
}
#endif /* SWIFT_TRACERS_EAGLE_IO_H */
