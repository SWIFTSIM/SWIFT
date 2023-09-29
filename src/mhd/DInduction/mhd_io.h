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
#ifndef SWIFT_DI_MHD_IO_H
#define SWIFT_DI_MHD_IO_H

#include "io_properties.h"
#include "statistics.h"

/**
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 * @return number of fields readed
 */
INLINE static int mhd_read_particles(struct part* parts,
                                     struct io_props* list) {

  list[0] = io_make_input_field("MagneticFluxDensity", FLOAT, 3, COMPULSORY,
                                UNIT_CONV_MAGNETIC_FIELD, parts,
                                mhd_data.BPred);  // CHECK XXX IF FULL STEP
  return 1;
}

INLINE static void convert_B(const struct engine* e, const struct part* p,
                             const struct xpart* xp, float* ret) {
  //  float a_fac = pow(e->cosmology->a, 3.f / 2.f * (hydro_gamma - 1.f) - 2.f);
  ret[0] =
      p->mhd_data.BPred[0];  // * sqrt(e->hydro_properties->mhd.mu0) * a_fac;
  ret[1] =
      p->mhd_data.BPred[1];  // * sqrt(e->hydro_properties->mhd.mu0) * a_fac;
  ret[2] =
      p->mhd_data.BPred[2];  // * sqrt(e->hydro_properties->mhd.mu0) * a_fac;
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

  list[0] = io_make_output_field("MagneticFluxDensity", FLOAT, 3,
                                 UNIT_CONV_MAGNETIC_FIELD, mhd_comoving_factor,
                                 parts, mhd_data.BPred,
                                 "Co-moving Magnetic field of the particles");

  list[1] = io_make_output_field(
      "MagneticDivergence", FLOAT, 1, UNIT_CONV_MAGNETIC_DIVERGENCE,
      mhd_comoving_factor - 1.f, parts, mhd_data.divB,
      "co-moving DivB of the particles");

  list[2] = io_make_output_field(
      "DednerScalar", FLOAT, 1, UNIT_CONV_ELECTRIC_CHARGE_FIELD_STRENGTH,
      mhd_comoving_factor + 1.f, parts, mhd_data.phi, "Dedner scalar associated to the particle");

  return 3;
}

/**
 * @brief Writes the current model of MHD to the file
 * @param h_grpsph The HDF5 group in which to write
 */
INLINE static void mhd_write_flavour(hid_t h_grpsph) {
  /* write XXX atributes for the implementation */
  /* really detail here */
  io_write_attribute_s(h_grpsph, "MHD Flavour",
                       "Federico - Direct Induction, divB substraction & "
                       "Dedner cleaning. Dolag & Stasyszyn (2009).");
}

#endif /* SWIFT_DI_MHD_IO_H */
