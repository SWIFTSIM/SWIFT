/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2026 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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
#ifndef SWIFT_STRENGTH_STATE_SETUP_H
#define SWIFT_STRENGTH_STATE_SETUP_H

/**
 * @file strength/strength_state_setup.h
 * @brief Setting up the strength state of a particle.
 */

#include "equation_of_state.h"

/**
 * @brief Sets the method-specific properties of the strength state from the
 * particle.
 *
 * @param s The strength state of the particle.
 * @param p The particle of interest.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void strength_state_set_method_fields(
    struct strength_state *s, const struct part *restrict p,
    const struct mat_params *mp) {

  yield_weakening_thermal_set_state(s, p, mp);
  yield_weakening_density_set_state(s, p, mp);
  yield_set_state(s, p, mp);
  damage_set_state(s, p, mp);
  artif_stress_set_state(s, p, mp);
}

/**
 * @brief Initialises the strength state of a particle.
 *
 * @param s The strength state to be initialised.
 * @param mat_id The material ID, for the EoS.
 * @param phase The phase.
 * @param mass The mass.
 * @param density The density.
 * @param u The specific internal energy.
 * @param damage The damage.
 * @param deviatoric_stress_tensor The deviatoric stress tensor.
 * @param strain_rate_tensor The strain rate tensor.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void strength_state_init(
    struct strength_state *s, const enum eos_planetary_material_id mat_id,
    const int phase, const float mass, const float density, const float u,
    const struct strength_damage damage,
    const struct sym_matrix deviatoric_stress_tensor,
    const struct sym_matrix strain_rate_tensor, const struct mat_params *mp) {

  /* Properties of the particle. */
  s->mat_id = mat_id;
  s->phase = phase;
  s->mass = mass;
  s->density = density;
  s->u = u;
  s->damage = damage;
  s->deviatoric_stress_tensor = deviatoric_stress_tensor;
  s->strain_rate_tensor = strain_rate_tensor;

  /* Quantities derived from the above are computed when first needed. */
  s->cache = (struct strength_state_cache){0};
}

/**
 * @brief Get the pressure, from the EoS.
 *
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float strength_state_get_pressure(
    struct strength_state *s, const struct mat_params *mp) {

  if (!s->cache.pressure_computed) {
    s->cache.pressure =
        gas_pressure_from_internal_energy(s->density, s->u, s->mat_id);
    s->cache.pressure_computed = true;
  }
  return s->cache.pressure;
}

/**
 * @brief Get the temperature, from the EoS.
 *
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float strength_state_get_temperature(
    struct strength_state *s, const struct mat_params *mp) {

  if (!s->cache.temperature_computed) {
    s->cache.temperature =
        gas_temperature_from_internal_energy(s->density, s->u, s->mat_id);
    s->cache.temperature_computed = true;
  }
  return s->cache.temperature;
}

/**
 * @brief Get the J_2 invariant of the deviatoric stress tensor.
 *
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float strength_state_get_J_2(
    struct strength_state *s, const struct mat_params *mp) {

  if (!s->cache.J_2_computed) {
    s->cache.J_2 =
        strength_compute_deviatoric_sym_matrix_J_2(s->deviatoric_stress_tensor);
    s->cache.J_2_computed = true;
  }
  return s->cache.J_2;
}

/**
 * @brief Get the stress tensor, including the effects of damage.
 *
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static struct sym_matrix
strength_state_get_stress_tensor(struct strength_state *s,
                                 const struct mat_params *mp) {

  if (!s->cache.stress_tensor_computed) {
    strength_compute_stress_tensor(&s->cache.stress_tensor, s, mp);
    s->cache.stress_tensor_computed = true;
  }
  return s->cache.stress_tensor;
}

#endif /* SWIFT_STRENGTH_STATE_SETUP_H */
