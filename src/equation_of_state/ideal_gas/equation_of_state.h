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

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "common_io.h"
#include "inline.h"
#include "physical_constants.h"

extern struct eos_parameters eos;

/**
 * @brief The parameters of the equation of state for the gas.
 *
 * This equation of state is parameter-free.
 */
struct eos_parameters {};

/**
 * @brief Returns the internal energy given density and entropy
 *
 * Computes \f$u = \frac{A\rho^{\gamma-1} }{\gamma - 1}\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$A\f$.
 */
__attribute__((always_inline, const)) INLINE static float
gas_internal_energy_from_entropy(float density, float entropy) {

  return entropy * pow_gamma_minus_one(density) *
         hydro_one_over_gamma_minus_one;
}

/**
 * @brief Returns the pressure given density and entropy
 *
 * Computes \f$P = A\rho^\gamma\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$A\f$.
 */
__attribute__((always_inline, const)) INLINE static float
gas_pressure_from_entropy(float density, float entropy) {

  return entropy * pow_gamma(density);
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
__attribute__((always_inline, const)) INLINE static float
gas_entropy_from_pressure(float density, float pressure) {

  return pressure * pow_minus_gamma(density);
}

/**
 * @brief Returns the sound speed given density and entropy
 *
 * Computes \f$c = \sqrt{\gamma A \rho^{\gamma-1}}\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$A\f$.
 */
__attribute__((always_inline, const)) INLINE static float
gas_soundspeed_from_entropy(float density, float entropy) {

  return sqrtf(hydro_gamma * pow_gamma_minus_one(density) * entropy);
}

/**
 * @brief Returns the entropy given density and internal energy
 *
 * Computes \f$A = \frac{(\gamma - 1)u}{\rho^{\gamma-1}}\f$.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline, const)) INLINE static float
gas_entropy_from_internal_energy(float density, float u) {

  return hydro_gamma_minus_one * u * pow_minus_gamma_minus_one(density);
}

/**
 * @brief Returns the pressure given density and internal energy
 *
 * Computes \f$P = (\gamma - 1)u\rho\f$.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline, const)) INLINE static float
gas_pressure_from_internal_energy(float density, float u) {

  return hydro_gamma_minus_one * u * density;
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
__attribute__((always_inline, const)) INLINE static float
gas_internal_energy_from_pressure(float density, float pressure) {

  return hydro_one_over_gamma_minus_one * pressure / density;
}

/**
 * @brief Returns the sound speed given density and internal energy
 *
 * Computes \f$c = \sqrt{\gamma (\gamma - 1) u }\f$.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline, const)) INLINE static float
gas_soundspeed_from_internal_energy(float density, float u) {

  return sqrtf(u * hydro_gamma * hydro_gamma_minus_one);
}

/**
 * @brief Returns the sound speed given density and pressure
 *
 * Computes \f$c = \sqrt{\frac{\gamma P}{\rho} }\f$.
 *
 * @param density The density \f$\rho\f$
 * @param P The pressure \f$P\f$
 */
__attribute__((always_inline, const)) INLINE static float
gas_soundspeed_from_pressure(float density, float P) {

  return sqrtf(hydro_gamma * P / density);
}

/**
 * @brief Initialize the eos parameters
 *
 * Nothing to do here since this EoS is parameter-free.
 *
 * @param e The #eos_parameters.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 */
INLINE static void eos_init(struct eos_parameters *e,
                            const struct phys_const *phys_const,
                            const struct unit_system *us,
                            struct swift_params *params) {}
/**
 * @brief Print the equation of state
 *
 * @param e The #eos_parameters
 */
INLINE static void eos_print(const struct eos_parameters *e) {

  message("Equation of state: Ideal gas.");

  message("Adiabatic index gamma: %f.", hydro_gamma);
}

#if defined(HAVE_HDF5)
/**
 * @brief Write equation of state information to the snapshot
 *
 * @param h_grpsph The HDF5 group in which to write
 * @param e The #eos_parameters
 */
INLINE static void eos_print_snapshot(hid_t h_grpsph,
                                      const struct eos_parameters *e) {

  io_write_attribute_f(h_grpsph, "Adiabatic index", hydro_gamma);

  io_write_attribute_s(h_grpsph, "Equation of state", "Ideal gas");
}
#endif

#endif /* SWIFT_IDEAL_GAS_EQUATION_OF_STATE_H */
