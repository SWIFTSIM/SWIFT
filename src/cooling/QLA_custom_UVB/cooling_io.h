/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *                    Stefan Arridge  (stefan.arridge@durham.ac.uk)
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
#ifndef SWIFT_COOLING_CONST_LAMBDA_IO_H
#define SWIFT_COOLING_CONST_LAMBDA_IO_H

/**
 * @file src/cooling/const_lambda/cooling_io.h
 * @brief i/o routines related to the "constant lambda" cooling function.
 *
 * This model assumes a constant cooling rate Lambda irrespective of redshift
 * or density.
 */

/* Config parameters. */
#include <config.h>

/* Local includes */
#include "cooling.h"
#include "engine.h"
#include "io_properties.h"

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of cooling to the file
 *
 * @param h_grp The HDF5 group in which to write
 * @param h_grp_columns The HDF5 group containing named columns
 * @param cooling The #cooling_function_data
 */
__attribute__((always_inline)) INLINE static void cooling_write_flavour(
    hid_t h_grp, hid_t h_grp_columns,
    const struct cooling_function_data* cooling) {

  io_write_attribute_s(h_grp, "Cooling Model", "Constant Lambda");
  io_write_attribute_d(h_grp, "Lambda/n_H^2 [cgs]", cooling->lambda_nH2_cgs);
  io_write_attribute_i(h_grp, "Rapid cooling", cooling->rapid_cooling);
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
 *
 * @return Returns the number of fields to write.
 */
__attribute__((always_inline)) INLINE static int cooling_write_particles(
    const struct part* parts, const struct xpart* xparts,
    struct io_props* list) {

  list[0] = io_make_output_field_convert_part(
      "Temperatures", FLOAT, 1, UNIT_CONV_TEMPERATURE, 0.f, parts, xparts,
      convert_part_T, "Temperatures of the gas particles");

  list[1] = io_make_output_field(
      "RadiatedEnergies", FLOAT, 1, UNIT_CONV_ENERGY, 0.f, xparts,
      cooling_data.radiated_energy,
      "Thermal energies radiated by the cooling mechanism");

  return 2;
}

#endif /* SWIFT_COOLING_CONST_LAMBDA_IO_H */
