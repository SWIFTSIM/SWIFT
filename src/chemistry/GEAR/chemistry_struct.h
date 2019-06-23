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
#ifndef SWIFT_CHEMISTRY_STRUCT_GEAR_H
#define SWIFT_CHEMISTRY_STRUCT_GEAR_H

/**
 * @brief The individual elements traced in the model.
 */
enum chemistry_element {
  chemistry_element_O = 0,
  chemistry_element_Mg,
  chemistry_element_S,
  chemistry_element_Fe,
  chemistry_element_Zn,
  chemistry_element_Sr,
  chemistry_element_Y,
  chemistry_element_Ba,
  chemistry_element_Eu,
  chemistry_element_count
};

/**
 * @brief Global chemical abundance information.
 */
struct chemistry_global_data {

  /* Initial metallicity Z */
  float initial_metallicity;
};

/**
 * @brief Properties of the chemistry function.
 */
struct chemistry_part_data {

  /*! Fraction of the particle mass in a given element */
  float metal_mass_fraction[chemistry_element_count];

  /*! Smoothed fraction of the particle mass in a given element */
  float smoothed_metal_mass_fraction[chemistry_element_count];

  float Z;
};

/**
 * @brief Chemical abundances traced by the #bpart in the GEAR model.
 */
struct chemistry_bpart_data {};

#endif /* SWIFT_CHEMISTRY_STRUCT_GEAR_H */
