/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Matthieu Schaller (schaller@strw.leideuniv.nl)
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
#ifndef SWIFT_DIRECT_INDUCTION_MHD_IO_H
#define SWIFT_DIRECT_INDUCTION_MHD_IO_H

#include "adiabatic_index.h"

/**
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 * @return number of fields readed
 */
INLINE static int mhd_read_particles(struct part* parts,
                                     struct io_props* list) {

  list[0] =
      io_make_input_field("MagneticFluxDensities", FLOAT, 3, COMPULSORY,
                          UNIT_CONV_MAGNETIC_FIELD, parts, mhd_data.B_over_rho);

  return 1;
}

INLINE static void convert_B(const struct engine* e, const struct part* p,
                             const struct xpart* xp, float* ret) {

  ret[0] = xp->mhd_data.B_over_rho_full[0] * p->rho;
  ret[1] = xp->mhd_data.B_over_rho_full[1] * p->rho;
  ret[2] = xp->mhd_data.B_over_rho_full[2] * p->rho;
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extended particle array.
 * @param list The list of i/o properties to write.
 * @return num_fields The number of i/o fields to write.
 */
INLINE static int mhd_write_particles(const struct part* parts,
                                      const struct xpart* xparts,
                                      struct io_props* list) {

  list[0] = io_make_output_field_convert_part(
      "MagneticFluxDensities", FLOAT, 3, UNIT_CONV_MAGNETIC_FIELD, -1.5f * hydro_gamma, parts,
      xparts, convert_B, "Magnetic flux densities of the particles");

  list[1] = io_make_output_field(
      "MagneticDivergences", FLOAT, 1, UNIT_CONV_MAGNETIC_DIVERGENCE, -1.5f * hydro_gamma - 1.f,
      parts, mhd_data.divB, "co-moving DivB  of the particle");

  list[2] = io_make_output_field(
      "DednerScalars", FLOAT, 1, UNIT_CONV_ELECTRIC_CHARGE_FIELD_STRENGTH, 1.5f * hydro_gamma + 1.f,
      parts, mhd_data.psi_over_ch, "Dedner scalar associated to the particle");

  list[3] = io_make_output_field(
      "DednerScalarsdt", FLOAT, 1,
      UNIT_CONV_ELECTRIC_CHARGE_FIELD_STRENGTH_RATE, 1.f, parts,
      mhd_data.psi_over_ch_dt,
      "Time derivative of Dedner scalar associated to particle");

  list[4] = io_make_output_field(
      "MagneticFluxDensitiesdt", FLOAT, 3, UNIT_CONV_MAGNETIC_FIELD_PER_TIME,
      1.f, parts, mhd_data.B_over_rho_dt,
      "Time derivative of Magnetic flux densities of the particles");

  list[5] = io_make_output_field(
      "MagneticFluxCurl", FLOAT, 3, UNIT_CONV_MAGNETIC_CURL, -1.5f * hydro_gamma - 1.f, parts,
      mhd_data.curl_B, "The curl of Magnetic flux densities of the particles");

  list[6] = io_make_output_field(
      "AlphaAR", FLOAT, 1, UNIT_CONV_NO_UNITS, 1.f, parts, mhd_data.alpha_AR,
      "Artificial resistivity switch of the particles");

  return 7;
}

/**
 * @brief Writes the current model of MHD to the file
 * @param h_grpsph The HDF5 group in which to write
 */
INLINE static void mhd_write_flavour(hid_t h_grpsph) {

  io_write_attribute_s(
      h_grpsph, "MHD Flavour",
      "Orestis - Direct Induction, divB Subtraction, "
      "Artificial Resistivity & Dedner Cleaning. Price et al. (2018).");
}

#endif /* SWIFT_DIRECT_INDUCTION_MHD_IO_H */
