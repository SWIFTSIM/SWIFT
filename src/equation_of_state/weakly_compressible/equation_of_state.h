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

//RGB this version is suitable for weakly compressible fluids typically used in engineering. This version is based on Morris, Fox & Zhu, 1997

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
 * RGB this version is suitable for weakly compressible fluids typically used in engineering. This version is based on Morris, Fox & Zhu, 1997
 */
struct eos_parameters {

  /*! Follows Morris et al Section 2.4, Eq 10. */
  float wc_soundspeed  /* this needs to be set by user to ensure density variations are sufficiently small */
};

/**
 * @brief Returns the internal energy given density and entropy
 *
 * INternal energy is not defined. Return an error. 
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float
gas_internal_energy_from_entropy(float density, float entropy) {
  error("internal energy is not defined in weakly copressible EOS");
  return 0.f;
}

/**
 * @brief Returns the pressure given density and entropy
 *
 * Since we are using an Morris et al EoS, the entropy value is ignored.
 * Computes \f$P = c^2 \rho\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float gas_pressure_from_entropy(
    float density, float entropy) {
  return eos.wc_soundspeed**2  * density ;
}

/**
 * @brief Returns the entropy given density and pressure.
 *
 * Since we are using an weakly compressible EoS, entropy is not defined. Return an error.
 *
 * @param density The density \f$\rho\f$.
 * @param pressure The pressure \f$P\f$ (ignored).
 * @return The entropy \f$A\f$.
 */
__attribute__((always_inline)) INLINE static float gas_entropy_from_pressure(
    float density, float pressure) {

  error("entropy is not defined in weakly copressible EOS");
  return 0.f;
}

/**
 * @brief Returns the sound speed given density and entropy
 *
 * Since we are using an weakly compressible EoS, this is just the parameter.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float gas_soundspeed_from_entropy(
    float density, float entropy) {
  return eos.wc_soundspeed ;
}

/**
 * @brief Returns the entropy given density and internal energy
 *
 *  Since we are using an weakly compressible EoS, entropy is not defined.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
gas_entropy_from_internal_energy(float density, float u) {

  // RGB entropy is not defined, so this is an error. 
  error("entropy is not defined in weakly copressible EOS");
  return 0.f ;
}

/**
 * @brief Returns the pressure given density and internal energy
 *
 * Since we are using an weakly compressible EOD, internal energy is ignored.
 * Computes \f$P = c^2 \rho\f$.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$ (ignored)
 */
__attribute__((always_inline)) INLINE static float
gas_pressure_from_internal_energy(float density, float u) {
  return eos.wc_soundspeed**2  * density ;
}

/**
 * @brief Returns the internal energy given density and pressure.
 *
 * Internal energy is not defined for the EOS. Returns an error.
 *
 * @param density The density \f$\rho\f$ (ignored).
 * @param pressure The pressure \f$P\f$ (ignored).
 * @return An error
 */
__attribute__((always_inline)) INLINE static float
gas_internal_energy_from_pressure(float density, float pressure) {
  error("internal energy is not defined in weakly copressible EOS");
  return 0.f ;
}

/**
 * @brief Returns the sound speed given density and internal energy
 *
 * Since we are using an weakly compressible EOS, this is just the parameter value.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$ (ignored)
 */
__attribute__((always_inline)) INLINE static float
gas_soundspeed_from_internal_energy(float density, float u) {
  return eos.wc_soundspeed;
}

/**
 * @brief Returns the sound speed given density and pressure
 *
 *  Since we are using an weakly compressible EOS, this is just the parameter value.
 *
 * @param density The density \f$\rho\f$ (ignored)
 * @param P The pressure \f$P\f$ (ignored)
 */
__attribute__((always_inline)) INLINE static float gas_soundspeed_from_pressure(
    float density, float P) {

  return  eos.wc_soundspeed;
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

  e->wc_soundspeed =
      parser_get_param_float(params, "EoS:wc_soundspeed");
}

/**
 * @brief Print the equation of state
 *
 * @param e The #eos_parameters
 */
__attribute__((always_inline)) INLINE static void eos_print(
    const struct eos_parameters *e) {
  message(
      "Equation of state: Weakly Compresssible EOS (Morris et al 1997). "
      "Numerical soundspeed per %f.",
      e->wc_soundspeed);
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
  io_write_attribute_s(h_grpsph, "equation of state", "Weakly Compressible, Morris et al 1997");
  io_write_attribute_f(h_grpsph, "Adiabatic index", hydro_gamma);
  io_write_attribute_f(h_grpsph, "numerical soundspeed",
                       e->wc_soundspeed);
}
#endif

#endif /* SWIFT_ISOTHERMAL_EQUATION_OF_STATE_H */
