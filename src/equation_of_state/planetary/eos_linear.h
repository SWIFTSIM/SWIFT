/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
 *               2024 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
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
#ifndef SWIFT_LINEAR_EQUATION_OF_STATE_H
#define SWIFT_LINEAR_EQUATION_OF_STATE_H

/**
 * @file equation_of_state/planetary/linear.h
 *
 * Contains the linear EOS functions for
 * equation_of_state/planetary/equation_of_state.h
 */

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "common_io.h"
#include "eos_setup.h"
#include "inline.h"
#include "physical_constants.h"
#include "utilities.h"

// Linear EoS parameters
struct linear_params {
  enum eos_planetary_material_id mat_id;
  float rho_0, c_s;
};

/*
    Read the parameters from a file.

    File contents
    -------------
    # header (2 lines)
    rho_0 (kg/m3)  c_s (m/s)
*/
INLINE static void set_linear_params(struct linear_params *linear,
                                     enum eos_planetary_material_id mat_id,
                                     char *param_file) {
  linear->mat_id = mat_id;

  // Load table contents from file
  FILE *f = fopen(param_file, "r");
  if (f == NULL) error("Failed to open the linear EoS file '%s'", param_file);

  // Skip header lines
  skip_lines(f, 2);

  // Read parameters (SI)
  int c;
  c = fscanf(f, "%f %f", &linear->rho_0, &linear->c_s);
  if (c != 2) error("Failed to read the linear EoS file %s", param_file);
}

// Convert to internal units
INLINE static void convert_units_linear(struct linear_params *linear,
                                        const struct unit_system *us) {

  struct unit_system si;
  units_init_si(&si);

  // SI to cgs
  linear->rho_0 *= units_cgs_conversion_factor(&si, UNIT_CONV_DENSITY);
  linear->c_s *= units_cgs_conversion_factor(&si, UNIT_CONV_SPEED);

  // cgs to internal
  linear->rho_0 /= units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  linear->c_s /= units_cgs_conversion_factor(us, UNIT_CONV_SPEED);
}

/**
 * @brief Returns the internal energy given density and entropy
 *
 * Computes \f$u = \frac{A\rho^{\gamma-1} }{\gamma - 1}\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$A\f$.
 */
INLINE static float linear_internal_energy_from_entropy(
    float density, float entropy, const struct linear_params *linear) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

/**
 * @brief Returns the pressure given density and entropy
 *
 * Computes \f$P = A\rho^\gamma\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$A\f$.
 */
INLINE static float linear_pressure_from_entropy(
    float density, float entropy, const struct linear_params *linear) {

  error("This EOS function is not yet implemented!");

  return 0.f;
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
INLINE static float linear_entropy_from_pressure(
    float density, float pressure, const struct linear_params *linear) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

/**
 * @brief Returns the sound speed given density and entropy
 *
 * Computes \f$c = \sqrt{\gamma A \rho^{\gamma-1}}\f$.
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$A\f$.
 */
INLINE static float linear_soundspeed_from_entropy(
    float density, float entropy, const struct linear_params *linear) {

  return linear->c_s;
}

/**
 * @brief Returns the entropy given density and internal energy
 *
 * Computes \f$A = (\gamma - 1) u \rho^{1-\gamma}\f$.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
INLINE static float linear_entropy_from_internal_energy(
    float density, float u, const struct linear_params *linear) {

#ifdef SWIFT_DEBUG_CHECKS
  warning("This EOS function is not yet implemented!");
#endif
  return 0.f;
}

/**
 * @brief Returns the pressure given density and internal energy
 *
 * Computes P = c_s^2 (rho - rho_0).
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
INLINE static float linear_pressure_from_internal_energy(
    float density, float u, const struct linear_params *linear) {

  return linear->c_s * linear->c_s * (density - linear->rho_0);
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
INLINE static float linear_internal_energy_from_pressure(
    float density, float pressure, const struct linear_params *linear) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

/**
 * @brief Returns the sound speed given density and internal energy
 *
 * Computes \f$c = \sqrt{\gamma (\gamma - 1) u }\f$.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
INLINE static float linear_soundspeed_from_internal_energy(
    float density, float u, const struct linear_params *linear) {

  return linear->c_s;
}

/**
 * @brief Returns the sound speed given density and pressure
 *
 * Computes \f$c = \sqrt{\frac{\gamma P}{\rho} }\f$.
 *
 * @param density The density \f$\rho\f$
 * @param P The pressure \f$P\f$
 */
INLINE static float linear_soundspeed_from_pressure(
    float density, float P, const struct linear_params *linear) {

  return linear->c_s;
}

// gas_temperature_from_internal_energy
INLINE static float linear_temperature_from_internal_energy(
    float density, float u, const struct linear_params *linear) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_density_from_pressure_and_temperature
INLINE static float linear_density_from_pressure_and_temperature(
    float P, float T, const struct linear_params *linear) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_density_from_pressure_and_internal_energy
INLINE static float linear_density_from_pressure_and_internal_energy(
    float P, float u, float rho_ref, float rho_sph,
    const struct linear_params *linear) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// material_phase_state_from_internal_energy
INLINE static float linear_phase_state_from_internal_energy(
    float density, float u, const struct mat_params *linear,
    const struct linear_params *linear_eos) {

  switch (linear->phase_state) {
    case mat_phase_state_fluid:
      return mat_phase_state_fluid;

    case mat_phase_state_solid:
      return mat_phase_state_solid;

    default:
      return mat_phase_state_fluid;
  }
}

#endif /* SWIFT_LINEAR_EQUATION_OF_STATE_H */
