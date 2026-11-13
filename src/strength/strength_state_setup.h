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

#include <float.h>

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

  yield_weakening_set_state(s, p, mp);
  yield_set_state(s, p, mp);
  damage_set_state(s, p, mp);
  artif_stress_set_state(s, p, mp);
}

/**
 * @brief Initialises the strength state of a particle.
 *
 * @param s The strength state to be initialised.
 * @param mat_id The material ID, for the EoS.
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
    const float mass, const float density, const float u,
    const struct strength_damage damage,
    const struct sym_matrix deviatoric_stress_tensor,
    const struct sym_matrix strain_rate_tensor, const struct mat_params *mp) {

  /* Properties of the particle. */
  s->mat_id = mat_id;
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
 * @brief Get the phase, from the material's state type and, for a variable
 * state type, its temperature and melting temperature.
 *
 * The choice between options is made here, so that each option can use
 * strength state quantities.
 *
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static enum mat_phase
strength_state_get_phase(struct strength_state *s, const struct mat_params *mp) {

  if (!s->cache.phase_computed) {
    switch (mp->state_type) {
      case mat_state_type_solid:
        s->cache.phase = mat_phase_solid;
        break;

      case mat_state_type_variable:
        s->cache.phase = (strength_state_get_temperature(s, mp) >
                          strength_state_get_T_melt(s, mp))
                             ? mat_phase_fluid
                             : mat_phase_solid;
        break;

      default:
        s->cache.phase = mat_phase_fluid;
    }
    s->cache.phase_computed = true;
  }
  return s->cache.phase;
}

/**
 * @brief Get the pressure, from the EoS, with tension limited by the minimum
 * pressure and by damage, depending on the yield stress method.
 *
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float strength_state_get_pressure(
    struct strength_state *s, const struct mat_params *mp) {

  if (!s->cache.pressure_computed) {
    float pressure =
        gas_pressure_from_internal_energy(s->density, s->u, s->mat_id);
    pressure = fmaxf(pressure, strength_state_get_P_min(s, mp));
    strength_modify_pressure_in_tension(&pressure, s, mp);
    s->cache.pressure = pressure;
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
 * @brief Get the melting temperature, from the material's scheme.
 *
 * The choice between options is made here, so that each option can use
 * strength state quantities.
 *
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float strength_state_get_T_melt(
    struct strength_state *s, const struct mat_params *mp) {

  if (!s->cache.T_melt_computed) {
    switch (mp->T_melt_scheme) {
      case T_melt_scheme_none:
        error("Material %d needs T_melt, but Material:T_melt_scheme is not set",
              s->mat_id);
        break;

      case T_melt_scheme_constant:
        s->cache.T_melt = material_T_melt_constant(mp);
        break;

      default:
        error("Unknown T_melt scheme %d", mp->T_melt_scheme);
    }
    s->cache.T_melt_computed = true;
  }
  return s->cache.T_melt;
}

/**
 * @brief Get the reference density, from the material's scheme.
 *
 * The choice between options is made here, so that each option can use
 * strength state quantities.
 *
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float strength_state_get_rho_0(
    struct strength_state *s, const struct mat_params *mp) {

  if (!s->cache.rho_0_computed) {
    switch (mp->rho_0_scheme) {
      case rho_0_scheme_none:
        error("Material %d needs rho_0, but Material:rho_0_scheme is not set",
              s->mat_id);
        break;

      case rho_0_scheme_constant:
        s->cache.rho_0 = material_rho_0_constant(mp);
        break;

      case rho_0_scheme_from_eos:
        s->cache.rho_0 = material_rho_0_from_eos(s->mat_id);
        break;

      default:
        error("Unknown rho_0 scheme %d", mp->rho_0_scheme);
    }
    s->cache.rho_0_computed = true;
  }
  return s->cache.rho_0;
}

/**
 * @brief Get the minimum pressure, from the material's scheme.
 *
 * The choice between options is made here, so that each option can use
 * strength state quantities.
 *
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float strength_state_get_P_min(
    struct strength_state *s, const struct mat_params *mp) {

  if (!s->cache.P_min_computed) {
    switch (mp->P_min_scheme) {
      case P_min_scheme_none:
        s->cache.P_min = -FLT_MAX;
        break;

      case P_min_scheme_constant:
        s->cache.P_min = material_P_min_constant(mp);
        break;

      case P_min_scheme_yield_envelope:
        s->cache.P_min = material_P_min_yield_envelope(mp, s->damage.total, s->mat_id);
        break;

      default:
        error("Unknown P_min scheme %d", mp->P_min_scheme);
    }
    s->cache.P_min_computed = true;
  }
  return s->cache.P_min;
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
    /* Damage scales the deviatoric stress tensor by the damage factor, and is
     * already applied to the pressure. */
    const float damage_factor = yield_damage_factor(s, mp);
    const float pressure = strength_state_get_pressure(s, mp);
    struct sym_matrix stress_tensor = s->deviatoric_stress_tensor;
    stress_tensor.xx = damage_factor * stress_tensor.xx - pressure;
    stress_tensor.yy = damage_factor * stress_tensor.yy - pressure;
    stress_tensor.zz = damage_factor * stress_tensor.zz - pressure;
    stress_tensor.xy *= damage_factor;
    stress_tensor.xz *= damage_factor;
    stress_tensor.yz *= damage_factor;
    s->cache.stress_tensor = stress_tensor;
    s->cache.stress_tensor_computed = true;
  }
  return s->cache.stress_tensor;
}

/**
 * @brief Get the principal stresses, the eigenvalues of the stress tensor.
 *
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 *
 * @return The three principal stresses.
 */
__attribute__((always_inline)) INLINE static const float *
strength_state_get_principal_stresses(struct strength_state *s,
                                      const struct mat_params *mp) {

  if (!s->cache.principal_stresses_computed) {
    sym_matrix_compute_eigenvalues(s->cache.principal_stresses,
                                   strength_state_get_stress_tensor(s, mp));
    s->cache.principal_stresses_computed = true;
  }
  return s->cache.principal_stresses;
}

/**
 * @brief Get the J_2 invariant of the strain rate tensor.
 *
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float strength_state_get_strain_rate_J_2(
    struct strength_state *s, const struct mat_params *mp) {

  if (!s->cache.strain_rate_J_2_computed) {
    s->cache.strain_rate_J_2 =
        strength_compute_sym_matrix_J_2(s->strain_rate_tensor);
    s->cache.strain_rate_J_2_computed = true;
  }
  return s->cache.strain_rate_J_2;
}

#ifdef STRENGTH_YIELD_STRESS
/**
 * @brief Get the yield stress, including damage and weakening.
 *
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float strength_state_get_yield_stress(
    struct strength_state *s, const struct mat_params *mp) {

  if (!s->cache.yield_stress_computed) {
    s->cache.yield_stress = yield_compute_yield_stress(s, mp);
    s->cache.yield_stress_computed = true;
  }
  return s->cache.yield_stress;
}
#endif

#endif /* SWIFT_STRENGTH_STATE_SETUP_H */
