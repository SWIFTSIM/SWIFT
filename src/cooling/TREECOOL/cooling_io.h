/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Matthieu Schaller (mschaller@lorentz.leidenuniv.nl)
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
#ifndef SWIFT_COOLING_TREECOOL_IO_H
#define SWIFT_COOLING_TREECOOL_IO_H

/**
 * @file src/cooling/TREECOOL/cooling_io.h
 * @brief i/o routines related to the TREECOOL cooling function (Katz,
 * Weinberg & Hernquist 1996).
 */

/* Config parameters. */
#include <config.h>

/* Local includes */
#include "cooling.h"
#include "engine.h"
#include "io_properties.h"

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of cooling to the file.
 *
 * @param h_grp The HDF5 group in which to write.
 * @param h_grp_columns The HDF5 group containing named columns.
 * @param cooling The #cooling_function_data.
 */
__attribute__((always_inline)) INLINE static void cooling_write_flavour(
    hid_t h_grp, hid_t h_grp_columns,
    const struct cooling_function_data *cooling) {

  io_write_attribute_s(h_grp, "Cooling Model", "TREECOOL");
  io_write_attribute_s(h_grp, "UV background table", cooling->TREECOOL_file);
  io_write_attribute_d(h_grp, "Hydrogen mass fraction", cooling->X_H);
  io_write_attribute_d(h_grp, "Helium mass fraction", cooling->Y_He);
  io_write_attribute_d(h_grp, "log10(T_min) [K]", cooling->log10_T_min);
  io_write_attribute_d(h_grp, "log10(T_max) [K]", cooling->log10_T_max);
  io_write_attribute_i(h_grp, "Compton cooling", cooling->with_Compton_cooling);
  io_write_attribute_f(h_grp, "Rapid cooling threshold",
                       cooling->rapid_cooling_threshold);
  io_write_attribute_f(h_grp, "UV background start redshift",
                       cooling->UV_background_start_redshift);
}
#endif

INLINE static void convert_part_T(const struct engine *e, const struct part *p,
                                  const struct xpart *xp, float *ret) {

  ret[0] = cooling_get_temperature(e->physical_constants, e->hydro_properties,
                                   e->internal_units, e->cosmology,
                                   e->cooling_func, p, xp);
}

/**
 * @brief Specifies which particle fields to write to a dataset.
 *
 * @param parts The particle array.
 * @param xparts The extended particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
__attribute__((always_inline)) INLINE static int cooling_write_particles(
    const struct part *parts, const struct xpart *xparts,
    struct io_props *list) {

  list[0] = io_make_output_field_convert_part(
      "Temperatures", FLOAT, 1, UNIT_CONV_TEMPERATURE, 0.f, parts, xparts,
      convert_part_T, "Temperatures of the gas particles");

  list[1] = io_make_physical_output_field(
      "RadiatedEnergies", FLOAT, 1, UNIT_CONV_ENERGY, 0.f, xparts,
      cooling_data.radiated_energy, /*can convert to comoving=*/0,
      "Thermal energies radiated by the cooling mechanism");

  list[2] = io_make_output_field(
      "ElectronFractions", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
      cooling_data.electron_fraction,
      "Electron number densities in units of the Hydrogen number densities, "
      "obtained by assuming ionization equilibrium");

  return 3;
}

#endif /* SWIFT_COOLING_TREECOOL_IO_H */
