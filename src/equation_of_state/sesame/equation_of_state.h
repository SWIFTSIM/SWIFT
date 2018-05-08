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
#ifndef SWIFT_SESAME_EQUATION_OF_STATE_H
#define SWIFT_SESAME_EQUATION_OF_STATE_H

/**
 * @file equation_of_state/sesame/equation_of_state.h
 *
 *          WIP PLACEHOLDER ONLY!
 */

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "common_io.h"
#include "inline.h"
#include "units.h"
#include "physical_constants.h"

extern struct eos_parameters eos;
/* ------------------------------------------------------------------------- */

// SESAME parameters
struct SESAME_params {
    int mat_id;
};

struct eos_parameters {
    struct SESAME_params SESAME_iron;
};

// Material identifier flags (material_ID = type_ID * type_factor + unit_ID)
#define type_factor 10
enum type_id {
    type_SESAME = 4
};
enum material_id {
    // SESAME
    SESAME_iron = type_SESAME*type_factor,
};

// Parameter values for each material (cgs units)
INLINE static void set_SESAME_iron(struct SESAME_params *mat) {
    mat->mat_id = SESAME_iron;
}

// Convert from cgs to internal units
INLINE static void convert_units_SESAME(
    struct SESAME_params *mat, const struct unit_system* us) {

}

/**
 * @brief Returns the internal energy given density and entropy
 *
 * NOT IMPLEMENTED!
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float
gas_internal_energy_from_entropy(float density, float entropy, int mat_id) {

  return 0;
}

/**
 * @brief Returns the pressure given density and entropy
 *
 * NOT IMPLEMENTED!
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float
gas_pressure_from_entropy(float density, float entropy, int mat_id) {

  return 0;
}

/**
 * @brief Returns the entropy given density and pressure.
 *
 * NOT IMPLEMENTED!
 *
 * @param density The density \f$\rho\f$.
 * @param pressure The pressure \f$P\f$.
 * @return The entropy \f$A\f$.
 */
__attribute__((always_inline)) INLINE static float
gas_entropy_from_pressure(float density, float pressure, int mat_id) {

  return 0;
}

/**
 * @brief Returns the sound speed given density and entropy
 *
 * NOT IMPLEMENTED!
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float
gas_soundspeed_from_entropy(float density, float entropy, int mat_id) {

  return 0;
}

/**
 * @brief Returns the entropy given density and internal energy
 *
 * NOT IMPLEMENTED!
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
gas_entropy_from_internal_energy(float density, float u, int mat_id) {

  return 0;
}

/**
 * @brief Returns the pressure given density and internal energy
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
gas_pressure_from_internal_energy(float density, float u, int mat_id) {

  return 0;
}

/**
 * @brief Returns the internal energy given density and pressure.
 *
 * NOT IMPLEMENTED!
 *
 * @param density The density \f$\rho\f$.
 * @param pressure The pressure \f$P\f$.
 * @return The internal energy \f$u\f$.
 */
__attribute__((always_inline)) INLINE static float
gas_internal_energy_from_pressure(float density, float pressure, int mat_id) {

  return 0;
}

/**
 * @brief Returns the sound speed given density and internal energy
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
gas_soundspeed_from_internal_energy(float density, float u, int mat_id) {

  return 0;
}

/**
 * @brief Returns the sound speed given density and pressure
 *
 * @param density The density \f$\rho\f$
 * @param P The pressure \f$P\f$
 */
__attribute__((always_inline)) INLINE static float
gas_soundspeed_from_pressure(float density, float P, int mat_id) {

  return 0;
}

/**
 * @brief Initialize the eos parameters
 *
 * @param e The #eos_parameters
 * @param params The parsed parameters
 */
__attribute__((always_inline)) INLINE static void eos_init(
    struct eos_parameters *e, const struct phys_const *phys_const,
    const struct unit_system *us, const struct swift_params *params) {

    // Set the parameters and load tables etc. for each material
    set_SESAME_iron(&e->SESAME_iron);

    // Convert to internal units
    convert_units_SESAME(&e->SESAME_iron, us);
}

/**
 * @brief Print the equation of state
 *
 * @param e The #eos_parameters
 */
__attribute__((always_inline)) INLINE static void eos_print(
    const struct eos_parameters *e) {

  message("Equation of state: SESAME.");
}

#if defined(HAVE_HDF5)
/**
 * @brief Write equation of state information to the snapshot
 *
 * @param h_grpsph The HDF5 group in which to write
 * @param e The #eos_parameters
 */
__attribute__((always_inline)) INLINE static void eos_print_snapshot(
    hid_t h_grpsph, const struct eos_parameters *e) {

  io_write_attribute_s(h_grpsph, "Equation of state", "SESAME");
}
#endif

#endif /* SWIFT_SESAME_EQUATION_OF_STATE_H */
