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
#ifndef SWIFT_CHEMISTRY_STRUCT_EAGLE_H
#define SWIFT_CHEMISTRY_STRUCT_EAGLE_H

/**
 * @brief The individual elements traced in the EAGLE model.
 */
enum eagle_chemisty_element {
  eagle_element_H = 0,
  eagle_element_He,
  eagle_element_C,
  eagle_element_N,
  eagle_element_O,
  eagle_element_Ne,
  eagle_element_Mg,
  eagle_element_Si,
  eagle_element_Fe,
  eagle_element_count
};

/**
 * @brief Chemical abundances traced in the EAGLE model.
 */
struct chemistry_part_data {

  /*! Fraction of the particle mass in a given element */
  float metal_mass_fraction[eagle_element_count];

  /*! Fraction of the particle mass in *all* metals */
  float metal_mass_fraction_total;

  /*! Smoothed fraction of the particle mass in a given element */
  float smoothed_metal_mass_fraction[eagle_element_count];

  /*! Smoothed fraction of the particle mass in *all* metals */
  float smoothed_metal_mass_fraction_total;

  /*! Mass coming from SNIa */
  float mass_from_SNIa;

  /*! Fraction of total gas mass in metals coming from SNIa */
  float metal_mass_fraction_from_SNIa;

  /*! Mass coming from AGB */
  float mass_from_AGB;

  /*! Fraction of total gas mass in metals coming from AGB */
  float metal_mass_fraction_from_AGB;

  /*! Mass coming from SNII */
  float mass_from_SNII;

  /*! Fraction of total gas mass in metals coming from SNII */
  float metal_mass_fraction_from_SNII;

  /*! Fraction of total gas mass in Iron coming from SNIa */
  float iron_mass_fraction_from_SNIa;

  /*! Smoothed fraction of total gas mass in Iron coming from SNIa */
  float smoothed_iron_mass_fraction_from_SNIa;
};

#endif /* SWIFT_CHEMISTRY_STRUCT_EAGLE_H */
