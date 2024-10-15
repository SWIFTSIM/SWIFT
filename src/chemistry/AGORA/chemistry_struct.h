/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023 Yves Revaz (yves.revaz@epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_STRUCT_AGORA_H
#define SWIFT_CHEMISTRY_STRUCT_AGORA_H

#define AGORA_CHEMISTRY_ELEMENT_COUNT \
  2  // we trace only two elements, the first one being null
#define AGORA_LABELS_SIZE 10

/**
 * @brief Global chemical abundance information.
 */
struct chemistry_global_data {

  /* Initial mass fraction */
  double initial_metallicities[AGORA_CHEMISTRY_ELEMENT_COUNT];

  /* Solar mass abundances read from the chemistry table */
  float solar_abundances[AGORA_CHEMISTRY_ELEMENT_COUNT];

  /*! Name of the different elements */
  char elements_name[AGORA_CHEMISTRY_ELEMENT_COUNT * AGORA_LABELS_SIZE];
};

/**
 * @brief Properties of the chemistry function for #part.
 */
struct chemistry_part_data {

  /*! Total mass of element in a particle. */
  double metal_mass[AGORA_CHEMISTRY_ELEMENT_COUNT];

#ifdef HYDRO_DOES_MASS_FLUX
  /*! Mass fluxes of the metals in a given element */
  double metal_mass_fluxes[AGORA_CHEMISTRY_ELEMENT_COUNT];
#endif

  /*! Smoothed fraction of the particle mass in a given element */
  double smoothed_metal_mass_fraction[AGORA_CHEMISTRY_ELEMENT_COUNT];
};

/**
 * @brief Properties of the chemistry function for #spart.
 */
struct chemistry_spart_data {

  /*! Fraction of the particle mass in a given element */
  double metal_mass_fraction[AGORA_CHEMISTRY_ELEMENT_COUNT];
};

/**
 * @brief Chemical abundances traced by the #bpart in the GEAR model.
 */
struct chemistry_bpart_data {};

/**
 * @brief Chemical abundances traced by the #sink in the GEAR model.
 */
struct chemistry_sink_data {};

#endif /* SWIFT_CHEMISTRY_STRUCT_AGORA_H */
