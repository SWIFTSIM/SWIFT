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
#ifndef SWIFT_ISOTHERMAL_EQUATION_OF_STATE_H
#define SWIFT_ISOTHERMAL_EQUATION_OF_STATE_H

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
 */
struct eos_parameters {

  /*! Thermal energy per unit mass */
  float isothermal_internal_energy;
};

/**
 * @brief Returns the internal energy given density and entropy
 *
 * Since we are using an isothermal EoS, the entropy and density values are
 * ignored.
 * Computes \f$u = u_{cst}\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float
gas_internal_energy_from_entropy(float density, float entropy) {

  return eos.isothermal_internal_energy;
}

/**
 * @brief Returns the pressure given density and entropy
 *
 * Since we are using an isothermal EoS, the entropy value is ignored.
 * Computes \f$P = (\gamma - 1)u_{cst}\rho\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float gas_pressure_from_entropy(
    float density, float entropy) {

  return hydro_gamma_minus_one * eos.isothermal_internal_energy * density;
}

/**
 * @brief Returns the entropy given density and pressure.
 *
 * Since we are using an isothermal EoS, the pressure value is ignored.
 * Computes \f$A = (\gamma - 1)u_{cst}\rho^{-(\gamma-1)}\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param pressure The pressure \f$P\f$ (ignored).
 * @return The entropy \f$A\f$.
 */
__attribute__((always_inline)) INLINE static float gas_entropy_from_pressure(
    float density, float pressure) {

  return hydro_gamma_minus_one * eos.isothermal_internal_energy *
         pow_minus_gamma_minus_one(density);
}

/**
 * @brief Returns the sound speed given density and entropy
 *
 * Since we are using an isothermal EoS, the entropy and density values are
 * ignored.
 * Computes \f$c = \sqrt{u_{cst} \gamma (\gamma-1)}\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float gas_soundspeed_from_entropy(
    float density, float entropy) {

  return sqrtf(eos.isothermal_internal_energy * hydro_gamma *
               hydro_gamma_minus_one);
}

/**
 * @brief Returns the entropy given density and internal energy
 *
 * Since we are using an isothermal EoS, the energy value is ignored.
 * Computes \f$S = \frac{(\gamma - 1)u_{cst}}{\rho^{\gamma-1}}\f$.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
gas_entropy_from_internal_energy(float density, float u) {

  return hydro_gamma_minus_one * eos.isothermal_internal_energy *
         pow_minus_gamma_minus_one(density);
}

/**
 * @brief Returns the pressure given density and internal energy
 *
 * Since we are using an isothermal EoS, the energy value is ignored.
 * Computes \f$P = (\gamma - 1)u_{cst}\rho\f$.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
gas_pressure_from_internal_energy(float density, float u) {

  return hydro_gamma_minus_one * eos.isothermal_internal_energy * density;
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
  return eos.isothermal_internal_energy;
}

/**
 * @brief Returns the sound speed given density and internal energy
 *
 * Since we are using an isothermal EoS, the energy and density values are
 * ignored.
 * Computes \f$c = \sqrt{u_{cst} \gamma (\gamma-1)}\f$.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
gas_soundspeed_from_internal_energy(float density, float u) {

  return sqrtf(eos.isothermal_internal_energy * hydro_gamma *
               hydro_gamma_minus_one);
}

/**
 * @brief Returns the sound speed given density and pressure
 *
 * Since we are using an isothermal EoS, the pressure and density values are
 * ignored.
 * Computes \f$c = \sqrt{u_{cst} \gamma (\gamma-1)}\f$.
 *
 * @param density The density \f$\rho\f$
 * @param P The pressure \f$P\f$
 */
__attribute__((always_inline)) INLINE static float gas_soundspeed_from_pressure(
    float density, float P) {

  return sqrtf(eos.isothermal_internal_energy * hydro_gamma *
               hydro_gamma_minus_one);
}

/**
 * @brief Initialize the eos parameters
 *
 * Read the constant internal energy from the parameter file.
 *
 * @param e The #eos_paramters.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 */
__attribute__((always_inline)) INLINE static void eos_init(
    struct eos_parameters *e, const struct phys_const *phys_const,
    const struct unit_system *us, struct swift_params *params) {

  e->isothermal_internal_energy =
      parser_get_param_float(params, "EoS:isothermal_internal_energy");
}

/**
 * @brief Print the equation of state
 *
 * @param e The #eos_parameters
 */
__attribute__((always_inline)) INLINE static void eos_print(
    const struct eos_parameters *e) {

  message(
      "Equation of state: Isothermal with internal energy "
      "per unit mass set to %f.",
      e->isothermal_internal_energy);

  message("Adiabatic index gamma: %f.", hydro_gamma);
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

  io_write_attribute_f(h_grpsph, "Adiabatic index", hydro_gamma);

  io_write_attribute_s(h_grpsph, "Equation of state", "Isothermal gas");
  io_write_attribute_f(h_grpsph, "Thermal energy per unit mass",
                       e->isothermal_internal_energy);
}
#endif

#endif /* SWIFT_ISOTHERMAL_EQUATION_OF_STATE_H */
