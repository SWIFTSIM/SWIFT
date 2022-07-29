/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Tsang Keung Chan (chantsangkeung@gmail.com)
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
#ifndef SWIFT_SPHM1RT_RT_SPECIES_AND_ELEMENTS_H
#define SWIFT_SPHM1RT_RT_SPECIES_AND_ELEMENTS_H

#include "inline.h"

/**
 * @file src/rt/SPHM1RT/rt_species_and_elements.h
 * @brief Main header file for species and elements definitions in SPHM1RT.
 */

/**
 * @brief The individual elements traced in the SPHM1RT model.
 */
enum rt_chemistry_element {
  rt_chemistry_element_H = 0,
  rt_chemistry_element_He,
  rt_chemistry_element_count
};

/**
 * @brief The individual species traced.
 */
enum rt_cooling_species {
  rt_sp_elec = 0, /* 0 */
  rt_sp_HI,       /* 1 */
  rt_sp_HII,      /* 2 */
  rt_sp_HeI,      /* 3 */
  rt_sp_HeII,     /* 4 */
  rt_sp_HeIII,    /* 5 */
  rt_species_count
};

/**
 * @brief Return a string containing the name of a given #rt_cooling_species.
 */
__attribute__((always_inline)) INLINE static const char* rt_get_species_name(
    enum rt_cooling_species spec) {

  static const char* rt_cooling_species_names[rt_species_count] = {
      "e", "HI", "HII", "HeI", "HeII", "HeIII"};

  return rt_cooling_species_names[spec];
}

/**
 * @brief Return a string containing the name of a given #rt_chemistry_element.
 */
__attribute__((always_inline)) INLINE static const char*
rt_chemistry_get_element_name(enum rt_chemistry_element elem) {

  static const char* rt_chemistry_element_names[rt_chemistry_element_count] = {
      "Hydrogen", "Helium"};

  return rt_chemistry_element_names[elem];
}

#endif /* SWIFT_SPHM1RT_RT_SPECIES_AND_ELEMENTS_H */
