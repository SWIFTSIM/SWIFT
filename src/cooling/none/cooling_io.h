/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_COOLING_NONE_IO_H
#define SWIFT_COOLING_NONE_IO_H

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "cooling.h"
#include "io_properties.h"

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of SPH to the file
 * @param h_grp The HDF5 group in which to write
 * @param cooling the parameters of the cooling function.
 */
__attribute__((always_inline)) INLINE static void cooling_write_flavour(
    hid_t h_grp, const struct cooling_function_data* cooling) {

  io_write_attribute_s(h_grp, "Cooling Model", "None");
}
#endif

INLINE static void convert_part_T(const struct engine* e, const struct part* p,
                                  const struct xpart* xp, float* ret) {

  ret[0] = cooling_get_temperature(e->physical_constants, e->hydro_properties,
                                   e->internal_units, e->cosmology,
                                   e->cooling_func, p, xp);
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extended particle array.
 * @param list The list of i/o properties to write.
 * @param cooling The #cooling_function_data
 *
 * @return Returns the number of fields to write.
 */
__attribute__((always_inline)) INLINE static int cooling_write_particles(
    const struct part* parts, const struct xpart* xparts, struct io_props* list,
    const struct cooling_function_data* cooling) {

  list[0] = io_make_output_field_convert_part("Temperature", FLOAT, 1,
                                              UNIT_CONV_TEMPERATURE, parts,
                                              xparts, convert_part_T);
  return 1;
}

#endif /* SWIFT_COOLING_NONE_IO_H */
