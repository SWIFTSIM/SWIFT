/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2025 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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
#ifndef SWIFT_DAMAGE_H
#define SWIFT_DAMAGE_H

/**
 * @file strength/damage/damage.h
 * @brief Core damage scheme that combines multiple fracture models.
 */

#include "const.h"
#include <float.h>
#include "math.h"

/**
 * @brief Get damage of particle at last drift time.
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float strength_get_damage(const struct part *restrict p) {

  return p->strength_data.damage;
}

/**
 * @brief Get damage of particle at last kick time.
 *
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static float strength_get_damage_full(const struct xpart *restrict xp) {

  return xp->strength_data.damage_full;
}

/**
 * @brief Set drift-time damage of particle.
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static void strength_set_damage(struct part *restrict p, const float damage) {

  p->strength_data.damage = damage;
}

/**
 * @brief Set kick-time damage of particle.
 *
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static void strength_set_damage_full(struct xpart *restrict xp, const float damage_full) {

  xp->strength_data.damage_full = damage_full;
}

/**
 * @brief Limits the time-step by the damage time-step.
 *
 * @param dt_cfl The hydro (+ strength) time-step.
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static void strength_compute_timestep_damage(
    float *dt_cfl, const struct part *restrict p) {

  *dt_cfl = fminf(*dt_cfl, p->strength_data.dt_damage);
}

/**
 * @brief Computes the damage-modified stress tensor.
 *
 * @param stress_tensor The stress tensor.
 * @param damaged_deviatoric_stress_tensor The damaged_deviatoric stress tensor, already modified by damage.
 * @param pressure The pressure.
 * @param damage The damage.
 */
__attribute__((always_inline)) INLINE static void damage_compute_stress_tensor(
    struct sym_matrix *stress_tensor, const struct sym_matrix damaged_deviatoric_stress_tensor, const float pressure, const float damage) {

  /* How damage affects the deviatoric stress tensor depends on the yield stress
   * method used. */
  *stress_tensor = damaged_deviatoric_stress_tensor;

  /* Damage weakens negative pressures so that fully damaged material cannot be
   * in tension. */
  const float effective_pressure = pressure > 0.f ? pressure : (1.f - damage) * pressure;
  stress_tensor->xx -= effective_pressure;
  stress_tensor->yy -= effective_pressure;
  stress_tensor->zz -= effective_pressure;
}

/**
 * @brief Sets the values of additional particle damage properties at a
 * kick time
 *
 * @param p The particle of interest.
 * @param xp The extended data of this particle.
 */
__attribute__((always_inline)) INLINE static void strength_reset_predicted_values_damage(
    struct part *restrict p, const struct xpart *restrict xp) {

  const float damage_full = strength_get_damage_full(xp);
  const float tensile_damage_full = damage_get_tensile_damage_full(xp);
  const float shear_damage_full = damage_get_shear_damage_full(xp);

  strength_set_damage(p, damage_full);
  damage_set_tensile_damage(p, tensile_damage_full);
  damage_set_shear_damage(p, shear_damage_full);
}

/**
 * @brief Evolves particle damage.
 *
 * @param damage The damage.
 * @param tensile_damage The tensile damage.
 * @param shear_damage The shear damage.
 * @param p The particle of interest.
 * @param stress_tensor The stress tensor.
 * @param strain_rate_tensor The strain rate tensor.
 * @param is_above_yield_criterion Whether the yield criterion is reached.
 * @param mat_id The material ID.
 * @param mass The particle mass.
 * @param density The density.
 * @param pressure The pressure.
 * @param dt_therm The time-step duration.
 */
__attribute__((always_inline)) INLINE static void damage_evolve(
    float *damage, float *tensile_damage, float *shear_damage, struct part *restrict p,
    const struct sym_matrix stress_tensor, const struct sym_matrix strain_rate_tensor,
    const int is_above_yield_criterion,
    const int mat_id, const float mass, const float density, const float pressure,
    const float dt_therm) {

  /* Evolve tensile damage. */
  damage_tensile_evolve(tensile_damage, p, stress_tensor, mat_id, mass, density, *damage, dt_therm);

  /* Evolve shear damage. */
  damage_shear_evolve(shear_damage, is_above_yield_criterion, strain_rate_tensor,
                      mat_id, pressure, dt_therm);

  /* Combine sources of damage. */
  *damage = fminf(*tensile_damage + *shear_damage, 1.f);
}

/**
 * @brief Evolves particle damage in the drift
 *
 * @param p The particle of interest.
 * @param stress_tensor The stress tensor.
 * @param strain_rate_tensor The strain rate tensor.
 * @param is_above_yield_criterion Whether the yield criterion is reached.
 * @param mat_id The material ID.
 * @param mass The particle mass.
 * @param density The density.
 * @param pressure The pressure.
 * @param dt_therm The time-step duration.
 */
__attribute__((always_inline)) INLINE static void damage_predict_evolve(
    struct part *restrict p, const struct sym_matrix stress_tensor,
    const struct sym_matrix strain_rate_tensor, const int is_above_yield_criterion,
    const int mat_id, const float mass, const float density, const float pressure,
    const float dt_therm) {

  /* Damage parameters set to values at drift time. */
  float damage = strength_get_damage(p);
  float tensile_damage = damage_get_tensile_damage(p);
  float shear_damage = damage_get_shear_damage(p);

  /* Evolve damage. */
  damage_evolve(&damage, &tensile_damage, &shear_damage, p, stress_tensor,
                strain_rate_tensor, is_above_yield_criterion,
                mat_id, mass, density, pressure, dt_therm);

  /* Update damage particle properties. */
  strength_set_damage(p, damage);
  damage_set_tensile_damage(p, tensile_damage);
  damage_set_shear_damage(p, shear_damage);
}

/**
 * @brief Evolves particle damage in the kick
 *
 * @param p The particle of interest.
 * @param xp The extended data of the particle of interest.
 * @param stress_tensor The stress tensor.
 * @param strain_rate_tensor The strain rate tensor.
 * @param is_above_yield_criterion Whether the yield criterion is reached.
 * @param mat_id The material ID.
 * @param mass The particle mass.
 * @param density The density.
 * @param pressure The pressure.
 * @param dt_therm The time-step duration.
 */
__attribute__((always_inline)) INLINE static void damage_kick_evolve(
    struct part *restrict p, struct xpart *restrict xp, const struct sym_matrix stress_tensor,
    const struct sym_matrix strain_rate_tensor, const int is_above_yield_criterion,
    const int mat_id, const float mass, const float density, const float pressure,
    const float dt_therm) {

  /* Damage parameters set to values at kick time. */
  float damage = strength_get_damage_full(xp);
  float tensile_damage = damage_get_tensile_damage_full(xp);
  float shear_damage = damage_get_shear_damage_full(xp);

  /* Evolve damage. */
  damage_evolve(&damage, &tensile_damage, &shear_damage, p, stress_tensor,
                strain_rate_tensor, is_above_yield_criterion,
                mat_id, mass, density, pressure, dt_therm);

  /* Update damage particle properties. */
  strength_set_damage_full(xp, damage);
  damage_set_tensile_damage_full(xp, tensile_damage);
  damage_set_shear_damage_full(xp, shear_damage);
}

/**
 * @brief Calculate the damage time-step.
 *
 * The time for damage to grow by damage_timestep_factor.
 *
 * @param p The particle of interest.
 * @param stress_tensor The stress tensor.
 * @param strain_rate_tensor The strain rate tensor.
 * @param is_above_yield_criterion Whether the yield criterion is reached.
 * @param mat_id The material ID.
 * @param mass The particle mass.
 * @param density The density.
 * @param pressure The pressure.
 */
__attribute__((always_inline)) INLINE static void damage_compute_timescale(
    struct part *restrict p, const struct sym_matrix stress_tensor,
    const struct sym_matrix strain_rate_tensor, const int is_above_yield_criterion,
    const int mat_id, const float mass, const float density, const float pressure) {

  /* Growth in damage used to set the time-step. */
  const float damage_timestep_factor = 0.1f; // ### Hardcoded for now. Treat this similarly to CFL in yaml

  /* Damage parameters set to values at drift time. */
  const float damage = strength_get_damage(p);
  const float tensile_damage = damage_get_tensile_damage(p);
  const float shear_damage = damage_get_shear_damage(p);

  float dt_damage = FLT_MAX;

  /* Further damage accumulation has no effect on fully damaged material. */
  if (damage >= 1.f) {
    p->strength_data.dt_damage = dt_damage;
    return;
  }

  /* Combined rate of total damage from both sources. */
  float dD_dt = 0.f;

  /* Tensile contribution. */
  if (tensile_damage < 1.f) {
    float tensile_cbrtD_dt = 0.f;
    int number_of_activated_flaws = 0;
    damage_tensile_compute_cbrtD_dt(&tensile_cbrtD_dt, &number_of_activated_flaws,
                                    p, stress_tensor, mat_id, mass, density,
                                    damage); // ### Should this be damage or tensile_damage?

    /* Tensile evolves cbrt(D), so convert its rate to a dD/dt via the growth in
     * cbrt(D) that corresponds to a growth of damage_timestep_factor in D. */
    if (tensile_cbrtD_dt > 0.f) {
      const float Delta_cbrtD =
          cbrtf(tensile_damage + damage_timestep_factor) - cbrtf(tensile_damage);
      dD_dt += damage_timestep_factor * tensile_cbrtD_dt / Delta_cbrtD;
    }
  }

  /* Shear contribution. The shear model evolves D directly. */
  float shear_dD_dt = 0.f;
  damage_shear_compute_dD_dt(&shear_dD_dt, is_above_yield_criterion,
                             strain_rate_tensor, mat_id, pressure, shear_damage);
  dD_dt += shear_dD_dt;

  /* Time for total damage to grow by damage_timestep_factor. */
  if (dD_dt > 0.f) {
    dt_damage = damage_timestep_factor / dD_dt;
  }

  p->strength_data.dt_damage = dt_damage;
}
/**
 * @brief Initialises the damage properties for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions or assignments between the particle
 * and extended particle fields.
 *
 * @param p The particle of interest.
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static void strength_first_init_part_damage(
    struct part *restrict p, struct xpart *restrict xp) {

  strength_set_damage(p, 0.f);
  damage_set_tensile_damage(p, 0.f);
  damage_set_shear_damage(p, 0.f);

  strength_set_damage_full(xp, 0.f);
  damage_set_tensile_damage_full(xp, 0.f);
  damage_set_shear_damage_full(xp, 0.f);
}

#endif /* SWIFT_DAMAGE_H */
