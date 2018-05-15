/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016   Matthieu Schaller (matthieu.schaller@durham.ac.uk).
 *               2018   Jacob Kegerreis (jacob.kegerreis@durham.ac.uk).
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
#ifndef SWIFT_ANEOS_EQUATION_OF_STATE_H
#define SWIFT_ANEOS_EQUATION_OF_STATE_H

/**
 * @file equation_of_state/planetary/aneos.h
 *
 * Contains the (M)ANEOS EOS functions for
 * equation_of_state/planetary/equation_of_state.h
 *
 * Adapted from the implementation in Gadget 2 of Cuk & Stewart (2012)
 *
 */

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "common_io.h"
#include "equation_of_state.h"
#include "inline.h"
#include "physical_constants.h"
#include "units.h"

// ANEOS parameters
struct ANEOS_params {
  enum eos_planetary_material_id mat_id;
};

// Parameter values for each material (cgs units)
INLINE static void set_ANEOS_iron(struct ANEOS_params *mat,
                                  enum eos_planetary_material_id mat_id) {
  mat->mat_id = mat_id;
}
INLINE static void set_MANEOS_forsterite(
    struct ANEOS_params *mat, enum eos_planetary_material_id mat_id) {
  mat->mat_id = mat_id;
}

// Convert from cgs to internal units
INLINE static void convert_units_ANEOS(struct ANEOS_params *mat,
                                       const struct unit_system *us) {}

// gas_internal_energy_from_entropy
INLINE static float ANEOS_internal_energy_from_entropy(
    float density, float entropy, const struct ANEOS_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

// gas_pressure_from_entropy
INLINE static float ANEOS_pressure_from_entropy(
    float density, float entropy, const struct ANEOS_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

// gas_entropy_from_pressure
INLINE static float ANEOS_entropy_from_pressure(
    float density, float pressure, const struct ANEOS_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

// gas_soundspeed_from_entropy
INLINE static float ANEOS_soundspeed_from_entropy(
    float density, float entropy, const struct ANEOS_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

// gas_entropy_from_internal_energy
INLINE static float ANEOS_entropy_from_internal_energy(
    float density, float u, const struct ANEOS_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

// gas_pressure_from_internal_energy
INLINE static float ANEOS_pressure_from_internal_energy(
    float density, float u, const struct ANEOS_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

// gas_internal_energy_from_pressure
INLINE static float ANEOS_internal_energy_from_pressure(
    float density, float P, const struct ANEOS_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

// gas_soundspeed_from_internal_energy
INLINE static float ANEOS_soundspeed_from_internal_energy(
    float density, float u, const struct ANEOS_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

// gas_soundspeed_from_pressure
INLINE static float ANEOS_soundspeed_from_pressure(
    float density, float P, const struct ANEOS_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

#endif /* SWIFT_ANEOS_EQUATION_OF_STATE_H */
