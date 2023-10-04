/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_STRUCT_GEAR_DIFFUSION_H
#define SWIFT_CHEMISTRY_STRUCT_GEAR_DIFFUSION_H

/**
 * @brief Global chemical abundance information.
 */
struct chemistry_global_data {

  /* Initial mass fraction */
  double initial_metallicities[GEAR_CHEMISTRY_ELEMENT_COUNT];

  /*! Diffusion coefficent (no unit) */
  float C;
};

/**
 * @brief Properties of the chemistry function for #part.
 */
struct chemistry_part_data {

  /*! Total mass of element in a particle.
    This field is available only outside the density hydro loop. */
  double metal_mass[GEAR_CHEMISTRY_ELEMENT_COUNT];

  /*! Smoothed fraction of the particle mass in a given element */
  double smoothed_metal_mass_fraction[GEAR_CHEMISTRY_ELEMENT_COUNT];

  /*! Diffusion coefficient */
  float diff_coef;

  /*! Variation of the metal mass */
  double metal_mass_dt[GEAR_CHEMISTRY_ELEMENT_COUNT];

  /*! shear tensor in internal and physical units. */
  float S[3][3];
};

/**
 * @brief Properties of the chemistry function for #spart.
 */
struct chemistry_spart_data {

  /*! Fraction of the particle mass in a given element */
  double metal_mass_fraction[GEAR_CHEMISTRY_ELEMENT_COUNT];
};

/**
 * @brief Chemical abundances traced by the #bpart in the GEAR model.
 */
struct chemistry_bpart_data {};

/**
 * @brief Chemical abundances traced by the #sink in the GEAR model.
 */
struct chemistry_sink_data {};

#endif /* SWIFT_CHEMISTRY_STRUCT_GEAR_DIFFUSION_H */
