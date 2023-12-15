/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023  Orestis Karapiperis (karapiperis@lorentz.leidenuniv.nl).
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
#ifndef SWIFT_BAROTROPIC_GAS_EQUATION_OF_STATE_H
#define SWIFT_BAROTROPIC_GAS_EQUATION_OF_STATE_H

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
 * Barotropic equation of state from Hennebelle et al., 2008, A&A, 477, 9
 * reimplemented following Pakmor et al., 2011, MNRAS, 418, 1392
 */
struct eos_parameters {

  /*! Square of barotropic sound speed in vacuum */
  float vacuum_sound_speed2;

  /*! Inverse of the core density */
  float inverse_core_density;
};

/**
 * @brief Returns the internal energy given density and entropy
 *
 * Since we are using a barotropic EoS, the entropy value is ignored.
 * Computes \f$u = c_0^2 \frac{1 + \left(\frac{\rho}{\rho_c}\right)^\gamma
 * }{\gamma - 1}\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$A\f$.
 */
__attribute__((always_inline, const)) INLINE static float
gas_internal_energy_from_entropy(float density, float entropy) {

  const float density_frac = density * eos.inverse_core_density;
  const float density_factor = pow_gamma(density_frac);

  return eos.vacuum_sound_speed2 * sqrtf(1.0f + density_factor) *
         hydro_one_over_gamma_minus_one;
}

/**
 * @brief Returns the pressure given density and entropy
 *
 * Since we are using a barotropic EoS, the entropy value is ignored.
 * Computes \f$P = c_0^2 \left(1 +
 * \left(\frac{\rho}{\rho_c}\right)^\gamma\right)\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$A\f$.
 */
__attribute__((always_inline, const)) INLINE static float
gas_pressure_from_entropy(float density, float entropy) {

  const float density_frac = density * eos.inverse_core_density;
  const float density_factor = pow_gamma(density_frac);

  return eos.vacuum_sound_speed2 * density * sqrtf(1.0f + density_factor);
}

/**
 * @brief Returns the entropy given density and pressure.
 *
 * Since we are using a barotropic EoS, the pressure value is ignored.
 * Computes \f$A = \rho^{1-\gamma}c_0^2 \left(1 +
 * \left(\frac{\rho}{\rho_c}\right)^\gamma\right)\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param pressure The pressure \f$P\f$.
 * @return The entropy \f$A\f$.
 */
__attribute__((always_inline, const)) INLINE static float
gas_entropy_from_pressure(float density, float pressure) {

  const float density_frac = density * eos.inverse_core_density;
  const float density_factor = pow_gamma(density_frac);

  return eos.vacuum_sound_speed2 * pow_minus_gamma_minus_one(density) *
         sqrtf(1.0f + density_factor);
}

/**
 * @brief Returns the sound speed given density and entropy
 *
 * Since we are using a barotropic EoS, the entropy is ignored.
 * Computes \f$c = \sqrt{c_0^2 \left(1 +
 * \left(\frac{\rho}{\rho_c}\right)^\gamma\right)}\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$A\f$.
 */
__attribute__((always_inline, const)) INLINE static float
gas_soundspeed_from_entropy(float density, float entropy) {

  const float density_frac = density * eos.inverse_core_density;
  const float density_factor = pow_gamma(density_frac);

  return sqrtf(eos.vacuum_sound_speed2 * sqrtf(1.0f + density_factor));
}

/**
 * @brief Returns the entropy given density and internal energy
 *
 * Since we are using a barotropic EoS, the internal energy value is ignored.
 * Computes \f$A = \rho^{1-\gamma}c_0^2 \left(1 +
 * \left(\frac{\rho}{\rho_c}\right)^\gamma\right)\f$.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline, const)) INLINE static float
gas_entropy_from_internal_energy(float density, float u) {

  const float density_frac = density * eos.inverse_core_density;
  const float density_factor = pow_gamma(density_frac);

  return eos.vacuum_sound_speed2 * pow_minus_gamma_minus_one(density) *
         sqrtf(1.0f + density_factor);
}

/**
 * @brief Returns the pressure given density and internal energy
 *
 * Since we are using a barotropic EoS, the internal energy value is ignored.
 * Computes \f$P = c_0^2 \left(1 +
 * \left(\frac{\rho}{\rho_c}\right)^\gamma\right)\f$.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline, const)) INLINE static float
gas_pressure_from_internal_energy(float density, float u) {

  const float density_frac = density * eos.inverse_core_density;
  const float density_factor = pow_gamma(density_frac);

  return eos.vacuum_sound_speed2 * density * sqrtf(1.0f + density_factor);
}

/**
 * @brief Returns the internal energy given density and pressure.
 *
 * Since we are using a barotropic EoS, the pressure value is ignored.
 * Computes \f$u = c_0^2 \frac{1 + \left(\frac{\rho}{\rho_c}\right)^\gamma
 * }{\gamma - 1}\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param pressure The pressure \f$P\f$.
 * @return The internal energy \f$u\f$.
 */
__attribute__((always_inline, const)) INLINE static float
gas_internal_energy_from_pressure(float density, float pressure) {

  const float density_frac = density * eos.inverse_core_density;
  const float density_factor = pow_gamma(density_frac);

  return eos.vacuum_sound_speed2 * sqrtf(1.0f + density_factor) *
         hydro_one_over_gamma_minus_one;
}

/**
 * @brief Returns the sound speed given density and internal energy
 *
 * Since we are using a barotropic EoS, the internal energy value is ignored.
 * Computes \f$c = \sqrt{c_0^2 \left(1 +
 * \left(\frac{\rho}{\rho_c}\right)^\gamma\right)}\f$.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline, const)) INLINE static float
gas_soundspeed_from_internal_energy(float density, float u) {

  const float density_frac = density * eos.inverse_core_density;
  const float density_factor = pow_gamma(density_frac);

  return sqrtf(eos.vacuum_sound_speed2 * sqrtf(1.0f + density_factor));
}

/**
 * @brief Returns the sound speed given density and pressure
 *
 * Since we are using a barotropic EoS, the pressure value is ignored.
 * Computes \f$c = \sqrt{c_0^2 \left(1 +
 * \left(\frac{\rho}{\rho_c}\right)^\gamma\right)}\f$.
 *
 * @param density The density \f$\rho\f$
 * @param P The pressure \f$P\f$
 */
__attribute__((always_inline, const)) INLINE static float
gas_soundspeed_from_pressure(float density, float P) {

  const float density_frac = density * eos.inverse_core_density;
  const float density_factor = pow_gamma(density_frac);

  return sqrtf(eos.vacuum_sound_speed2 * sqrtf(1.0f + density_factor));
}

/**
 * @brief Initialize the eos parameters
 *
 * Read the vacuum sound speed and core density from the parameter file.
 *
 * @param e The #eos_parameters.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 */
INLINE static void eos_init(struct eos_parameters *e,
                            const struct phys_const *phys_const,
                            const struct unit_system *us,
                            struct swift_params *params) {

  const float vacuum_sound_speed =
      parser_get_param_float(params, "EoS:barotropic_vacuum_sound_speed");
  e->vacuum_sound_speed2 = vacuum_sound_speed * vacuum_sound_speed;
  e->inverse_core_density =
      1. / parser_get_param_float(params, "EoS:barotropic_core_density");
}
/**
 * @brief Print the equation of state
 *
 * @param e The #eos_parameters
 */
INLINE static void eos_print(const struct eos_parameters *e) {

  message(
      "Equation of state: Barotropic gas with vacuum sound speed set to %f and "
      "core density set to %f.",
      sqrtf(e->vacuum_sound_speed2), 1. / e->inverse_core_density);

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

  io_write_attribute_s(h_grpsph, "Equation of state", "Barotropic gas");
}
#endif

#endif /* SWIFT_BAROTROPIC_GAS_EQUATION_OF_STATE_H */
