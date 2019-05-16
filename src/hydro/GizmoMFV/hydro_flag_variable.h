/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

#ifndef SWIFT_HYDRO_GIZMO_MFV_FLAG_VARIABLE_H
#define SWIFT_HYDRO_GIZMO_MFV_FLAG_VARIABLE_H

#include "io_properties.h"

enum gizmo_flag_variables {
  GIZMO_FLAG_ENTROPY = 0,
  GIZMO_FLAG_MORE_NEIGHBOURS = 1,
  GIZMO_FLAG_REVERT_TO_SPH = 2,
  GIZMO_FLAG_EKIN_SWITCH = 3,
  GIZMO_FLAG_GRAVITY_SWITCH = 4
};

#if defined(GIZMO_FLAG_VARIABLE)

/**
 * @brief Initialize the flag variable.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_flag_variable_init(
    struct part* p) {
  p->flag_variable = 0;
}

/**
 * @brief Add the flag variable to the snapshot output.
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 */
__attribute__((always_inline)) INLINE static void hydro_write_flag_variable(
    const struct part* parts, struct io_props* list, int* num_fields) {

  list[*num_fields] = io_make_output_field(
      "FlagVariable", INT, 1, UNIT_CONV_NO_UNITS, parts, flag_variable);
  ++(*num_fields);
}

/**
 * @brief Add the given flag value to the flag variable.
 *
 * @param p Particle.
 * @param value Flag value to add.
 */
__attribute__((always_inline)) INLINE static void hydro_add_flag(struct part* p,
                                                                 int value) {

  p->flag_variable |= (1 << value);
}

#else

/**
 * @brief Initialize the flag variable.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_flag_variable_init(
    struct part* p) {}

/**
 * @brief Add the flag variable to the snapshot output.
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 */
__attribute__((always_inline)) INLINE static void hydro_write_flag_variable(
    const struct part* parts, struct io_props* list, int* num_fields) {}

/**
 * @brief Add the given flag value to the flag variable.
 *
 * @param p Particle.
 * @param value Flag value to add.
 */
__attribute__((always_inline)) INLINE static void hydro_add_flag(struct part* p,
                                                                 int value) {}

#endif

#endif /* SWIFT_HYDRO_GIZMO_MFV_FLAG_VARIABLE_H */
