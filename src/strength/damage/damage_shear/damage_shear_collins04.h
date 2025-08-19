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
#ifndef SWIFT_DAMAGE_SHEAR_COLLINS04_H
#define SWIFT_DAMAGE_SHEAR_COLLINS04_H

/**
 * @file strength/damage/damage_shear/damage_shear_collins04.h
 * Collins+2004 shear fracture model.
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Get shear damage of particle at last drift time.
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float damage_get_shear_damage(const struct part *restrict p) {

  return p->strength_data.shear_damage;
}

/**
 * @brief Get shear damage of particle at last kick time.
 *
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static float damage_get_shear_damage_full(const struct xpart *restrict xp) {

  return xp->strength_data.shear_damage_full;
}

/**
 * @brief Set drift-time shear damage of particle.
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static void damage_set_shear_damage(struct part *restrict p, const float shear_damage) {

  p->strength_data.shear_damage = shear_damage;
}

/**
 * @brief Set kick-time shear damage of particle.
 *
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static void damage_set_shear_damage_full(struct xpart *restrict xp, const float shear_damage_full) {

  xp->strength_data.shear_damage_full = shear_damage_full;
}

/**
 * @brief Calculates the rate of damage accumulation due to shear
 *
 * Method parameters needed in material parameter file:
 * DamageShearCollins:
 *     brittle_to_ductile_pressure: brittle--ductile transition pressure (Pa).
 *     brittle_to_plastic_pressure: brittle--plastic transition pressure (Pa).
 *
 * @param shear_dD_dt The rate of shear damage accumulation.
 * @param deviatoric_stress_tensor The stress deviatoric tensor.
 * @param mat_id The material ID.
 * @param density The density.
 * @param u The specific internal energy.
 * @param yield_stress The yield stress.
 * @param shear_damage The shear damage.
 */
__attribute__((always_inline)) INLINE static void damage_shear_compute_dD_dt(
    float *shear_dD_dt, const struct sym_matrix deviatoric_stress_tensor, const int mat_id, const float density, const float u, const float yield_stress, const float shear_damage) {

  *shear_dD_dt = 0.f;

  /* Return dD/dt = 0 if already at max shear damage. */
  if (shear_damage >= 1.f) {
    return;
  }

  // ### I think this is wrong and I need to calculate a quantity that gets
  // ### evolved in time based on adding an invariant of the *plastic strain*
  // ### each timestep. Instead of current invariant of stress. I need to evolve
  // ### stain epsilon = strain_rate * dt in time, limited by Y/sqrtf(J_2) and
  // ### then the plastic strain is the accumulation of the strain that gets
  // ### reset based on Y/sqrtf(J_2). In practice this becomes:

  // ### Evolve *strain* based on strain rate.
  // ### Make sure to do evolution in co-rotating frame.
  // ### deviatoric strain gets reduced by yield criterion in the same way as deviatoric stress.
  // ### The accumulation of the amount that the deviatoric strain gets reduced by in each elements gets stored as the plastic strain.
  // ### The equivalent of the J_2 invariant of this plastic strain is added up each timestep to get a quantity epsilon_tot (Collins Eqn A6).

  /* Calculate pressure and J_2 invariant of stress tensor. */
  const float pressure =
        gas_pressure_from_internal_energy(density, u, mat_id);
  const float J_2 = strength_compute_stress_tensor_J_2(deviatoric_stress_tensor);
  const float strain_rate_invariant = sqrtf(J_2);

  /* Method parameters. */
  const float brittle_to_ductile_pressure = material_brittle_to_ductile_pressure(mat_id);
  const float brittle_to_plastic_pressure = material_brittle_to_plastic_pressure(mat_id);

  /* Compute plastic strain at failure. */
  float slope, intercept;
  if (pressure <= 0.f) {
    // ### In paper its unclear whether shear damage can still accumulate when
    // ### p < 0. However, in the iSALE manual Eqn. 4.81 it looks like this is right?
    slope = 0.f;
    intercept = 0.01f;

  } else if (pressure < brittle_to_ductile_pressure) {
    /* Transitions between 0.01 and 0.05 in brittle regime. */
    slope = (0.05f - 0.01f) / brittle_to_ductile_pressure;
    intercept = 0.01f;

  } else if (pressure < brittle_to_plastic_pressure) {
    /* Transitions between 0.05 and 0.1 in semi-brittle regime. */
    slope = (0.1f - 0.05f) / (brittle_to_plastic_pressure - brittle_to_ductile_pressure);
    intercept = 0.05f - slope * brittle_to_ductile_pressure;

  } else {
    /* Transitions between 0.1 and 1 in plastic regime. */
    // ### In paper its unclear whether it continues to increase past 1, like I
    // ### have here. But in the iSALE amunal Eqn. 4.83 it looks like it does.
    // ### Although with a different equation??
    slope = (1.f - 0.1f) / brittle_to_plastic_pressure;
    intercept = 0.1f - slope * brittle_to_plastic_pressure;
  }

  const float plastic_strain_at_failure = slope * pressure + intercept;

  // ### Main questions here are:
  // ### 1) what the eqn is in the "else" above:
  // ### 2) Can damage decrease if e.g.  plastic_strain_at_failure increases but evolved plastic strain is const

  // ### This is wrong. I need to calcualte based on the new evolved strain rate invariant
  *shear_dD_dt = (strain_rate_invariant / plastic_strain_at_failure);
}

/**
 * @brief Steps particle shear damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void damage_shear_apply_timestep_to_shear_damage(
    float *shear_damage, const float shear_dD_dt, const float dt_therm) {

  *shear_damage += fminf(shear_dD_dt * dt_therm, 1.f);
}

/**
 * @brief Evolves particle shear damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void damage_shear_evolve(
    float *shear_damage, struct part *restrict p, const struct sym_matrix deviatoric_stress_tensor,
    const int mat_id, const float density, const float u, const float yield_stress, const float dt_therm) {

  float shear_dD_dt;
  damage_shear_compute_dD_dt(&shear_dD_dt, deviatoric_stress_tensor, mat_id, density, u, yield_stress, *shear_damage);
  damage_shear_apply_timestep_to_shear_damage(shear_damage, shear_dD_dt, dt_therm);
}

#endif /* SWIFT_DAMAGE_SHEAR_COLLINS04_H */
