/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016   Matthieu Schaller (schaller@strw.leidenuniv.nl).
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
#ifndef SWIFT_IDEAL_GAS_EQUATION_OF_STATE_H
#define SWIFT_IDEAL_GAS_EQUATION_OF_STATE_H

/**
 * @file equation_of_state/planetary/ideal_gas.h
 *
 * Contains the ideal gas EOS functions for
 * equation_of_state/planetary/equation_of_state.h
 */

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "common_io.h"
#include "inline.h"
#include "physical_constants.h"

// Ideal gas parameters
struct idg_params {
  enum eos_planetary_material_id mat_id;
  float gamma, one_over_gamma_minus_one;
};

// Parameter values for each material
INLINE static void set_idg_def(struct idg_params *mat,
                               enum eos_planetary_material_id mat_id) {
  mat->mat_id = mat_id;
  mat->gamma = hydro_gamma;  // set by --with-adiabatic-index, default 5/3
  mat->one_over_gamma_minus_one = 1.f / (mat->gamma - 1.f);
}

/**
 * @brief Returns the internal energy given density and entropy
 *
 * Computes \f$u = \frac{A\rho^{\gamma-1} }{\gamma - 1}\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$A\f$.
 */
INLINE static float idg_internal_energy_from_entropy(
    float density, float entropy, const struct idg_params *mat) {

  return entropy * powf(density, mat->gamma - 1.f) *
         mat->one_over_gamma_minus_one;
}

/**
 * @brief Returns the pressure given density and entropy
 *
 * Computes \f$P = A\rho^\gamma\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$A\f$.
 */
INLINE static float idg_pressure_from_entropy(float density, float entropy,
                                              const struct idg_params *mat) {

  return entropy * powf(density, mat->gamma);
}

/**
 * @brief Returns the entropy given density and pressure.
 *
 * Computes \f$A = \frac{P}{\rho^-\gamma}\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param pressure The pressure \f$P\f$.
 * @return The entropy \f$A\f$.
 */
INLINE static float idg_entropy_from_pressure(float density, float pressure,
                                              const struct idg_params *mat) {

  return pressure * powf(density, -mat->gamma);
}

/**
 * @brief Returns the sound speed given density and entropy
 *
 * Computes \f$c = \sqrt{\gamma A \rho^{\gamma-1}}\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$A\f$.
 */
INLINE static float idg_soundspeed_from_entropy(float density, float entropy,
                                                const struct idg_params *mat) {

  return sqrtf(mat->gamma * powf(density, mat->gamma - 1.f) * entropy);
}

/**
 * @brief Returns the entropy given density and internal energy
 *
 * Computes \f$A = (\gamma - 1) u \rho^{1-\gamma}\f$.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
INLINE static float idg_entropy_from_internal_energy(
    float density, float u, const struct idg_params *mat) {

  return (mat->gamma - 1.f) * u * powf(density, 1.f - mat->gamma);
}

/**
 * @brief Returns the pressure given density and internal energy
 *
 * Computes \f$P = (\gamma - 1)u\rho\f$.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
INLINE static float idg_pressure_from_internal_energy(
    float density, float u, const struct idg_params *mat) {

  return (mat->gamma - 1.f) * u * density;
}

/**
 * @brief Returns the internal energy given density and pressure.
 *
 * Computes \f$u = \frac{1}{\gamma - 1}\frac{P}{\rho}\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param pressure The pressure \f$P\f$.
 * @return The internal energy \f$u\f$.
 */
INLINE static float idg_internal_energy_from_pressure(
    float density, float pressure, const struct idg_params *mat) {

  return mat->one_over_gamma_minus_one * pressure / density;
}

/**
 * @brief Returns the sound speed given density and internal energy
 *
 * Computes \f$c = \sqrt{\gamma (\gamma - 1) u }\f$.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
INLINE static float idg_soundspeed_from_internal_energy(
    float density, float u, const struct idg_params *mat) {

  return sqrtf(u * mat->gamma * (mat->gamma - 1.f));
}

/**
 * @brief Returns the sound speed given density and pressure
 *
 * Computes \f$c = \sqrt{\frac{\gamma P}{\rho} }\f$.
 *
 * @param density The density \f$\rho\f$
 * @param P The pressure \f$P\f$
 */
INLINE static float idg_soundspeed_from_pressure(float density, float P,
                                                 const struct idg_params *mat) {

  return sqrtf(mat->gamma * P / density);
}

#endif /* SWIFT_IDEAL_GAS_EQUATION_OF_STATE_H */
