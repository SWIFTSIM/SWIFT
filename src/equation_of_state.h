/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016   Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#ifndef SWIFT_EQUATION_OF_STATE_H
#define SWIFT_EQUATION_OF_STATE_H

/**
 * @file equation_of_state.h
 * @brief Defines the equation of state of the gas we simulate in the form of
 * relations between thermodynamic quantities. These are later used internally
 * by all hydro schemes
 */

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "debug.h"
#include "inline.h"

/* ------------------------------------------------------------------------- */
#if defined(EOS_IDEAL_GAS)

/**
 * @brief Returns the internal energy given density and entropy
 *
 * Computes \f$u = \frac{S\rho^{\gamma-1} }{\gamma - 1}\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float
gas_internal_energy_from_entropy(float density, float entropy) {

  return entropy * pow_gamma_minus_one(density) *
         hydro_one_over_gamma_minus_one;
}

/**
 * @brief Returns the pressure given density and entropy
 *
 * Computes \f$P = S\rho^\gamma\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float gas_pressure_from_entropy(
    float density, float entropy) {

  return entropy * pow_gamma(density);
}

/**
 * @brief Returns the entropy given density and pressure.
 *
 * Computes \f$A = \frac{P}{\rho^\gamma}\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param pressure The pressure \f$P\f$.
 * @return The entropy \f$A\f$.
 */
__attribute__((always_inline)) INLINE static float gas_entropy_from_pressure(
    float density, float pressure) {

  return pressure * pow_minus_gamma(density);
}

/**
 * @brief Returns the sound speed given density and entropy
 *
 * Computes \f$c = \sqrt{\gamma S \rho^{\gamma-1}}\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float gas_soundspeed_from_entropy(
    float density, float entropy) {

  return sqrtf(hydro_gamma * pow_gamma_minus_one(density) * entropy);
}

/**
 * @brief Returns the entropy given density and internal energy
 *
 * Computes \f$S = \frac{(\gamma - 1)u}{\rho^{\gamma-1}}\f$.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
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
__attribute__((always_inline)) INLINE static float
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
__attribute__((always_inline)) INLINE static float
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
__attribute__((always_inline)) INLINE static float
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
__attribute__((always_inline)) INLINE static float gas_soundspeed_from_pressure(
    float density, float P) {

  const float density_inv = 1.f / density;
  return sqrtf(hydro_gamma * P * density_inv);
}

/* ------------------------------------------------------------------------- */
#elif defined(EOS_ISOTHERMAL_GAS)

/**
 * @brief Returns the internal energy given density and entropy
 *
 * Since we are using an isothermal EoS, the entropy value is ignored
 * Computes \f$u = u_{cst}\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float
gas_internal_energy_from_entropy(float density, float entropy) {

  return const_isothermal_internal_energy;
}
/**
 * @brief Returns the pressure given density and entropy
 *
 * Since we are using an isothermal EoS, the entropy value is ignored
 * Computes \f$P = (\gamma - 1)u_{cst}\rho\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float gas_pressure_from_entropy(
    float density, float entropy) {

  return hydro_gamma_minus_one * const_isothermal_internal_energy * density;
}

/**
 * @brief Returns the entropy given density and pressure.
 *
 * Computes \f$A = \frac{P}{\rho^\gamma}\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param pressure The pressure \f$P\f$ (ignored).
 * @return The entropy \f$A\f$.
 */
__attribute__((always_inline)) INLINE static float gas_entropy_from_pressure(
    float density, float pressure) {

  return hydro_gamma_minus_one * const_isothermal_internal_energy *
         pow_minus_gamma_minus_one(density);
}

/**
 * @brief Returns the sound speed given density and entropy
 *
 * Since we are using an isothermal EoS, the entropy value is ignored
 * Computes \f$c = \sqrt{u_{cst} \gamma \rho^{\gamma-1}}\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float gas_soundspeed_from_entropy(
    float density, float entropy) {

  return sqrtf(const_isothermal_internal_energy * hydro_gamma *
               hydro_gamma_minus_one);
}

/**
 * @brief Returns the entropy given density and internal energy
 *
 * Since we are using an isothermal EoS, the energy value is ignored
 * Computes \f$S = \frac{(\gamma - 1)u_{cst}}{\rho^{\gamma-1}}\f$.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
gas_entropy_from_internal_energy(float density, float u) {

  return hydro_gamma_minus_one * const_isothermal_internal_energy *
         pow_minus_gamma_minus_one(density);
}

/**
 * @brief Returns the pressure given density and internal energy
 *
 * Since we are using an isothermal EoS, the energy value is ignored
 * Computes \f$P = (\gamma - 1)u_{cst}\rho\f$.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
gas_pressure_from_internal_energy(float density, float u) {

  return hydro_gamma_minus_one * const_isothermal_internal_energy * density;
}

/**
 * @brief Returns the internal energy given density and pressure.
 *
 * Just returns the constant internal energy.
 *
 * @param density The density \f$\rho\f$ (ignored).
 * @param pressure The pressure \f$P\f$ (ignored).
 * @return The internal energy \f$u\f$ (which is constant).
 */
__attribute__((always_inline)) INLINE static float
gas_internal_energy_from_pressure(float density, float pressure) {
  return const_isothermal_internal_energy;
}

/**
 * @brief Returns the sound speed given density and internal energy
 *
 * Since we are using an isothermal EoS, the energy value is ignored
 * Computes \f$c = \sqrt{u_{cst} \gamma \rho^{\gamma-1}}\f$.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
gas_soundspeed_from_internal_energy(float density, float u) {

  return sqrtf(const_isothermal_internal_energy * hydro_gamma *
               hydro_gamma_minus_one);
}

/**
 * @brief Returns the sound speed given density and pressure
 *
 * Since we are using an isothermal EoS, the pressure value is ignored
 * Computes \f$c = \sqrt{u_{cst} \gamma \rho^{\gamma-1}}\f$.
 *
 * @param density The density \f$\rho\f$
 * @param P The pressure \f$P\f$
 */
__attribute__((always_inline)) INLINE static float gas_soundspeed_from_pressure(
    float density, float P) {

  return sqrtf(const_isothermal_internal_energy * hydro_gamma *
               hydro_gamma_minus_one);
}

/* ------------------------------------------------------------------------- */
#else

#error "An Equation of state needs to be chosen in const.h !"

#endif

#endif /* SWIFT_EQUATION_OF_STATE_H */
