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
 * @brief Calculates the damage accumulated due to shear
 *
 * Method parameters needed in material parameter file:
 * DamageShearCollins:
 *     brittle_to_ductile_pressure: brittle--ductile transition pressure (Pa).
 *     brittle_to_plastic_pressure: brittle--plastic transition pressure (Pa).
 *
 * @param shear_damage The updated shear damage.
 * @param p The particle of interest.
 * @param mat_id The material ID.
 * @param density The density.
 * @param u The specific internal energy.
 */
__attribute__((always_inline)) INLINE static void damage_shear_evolve(
    float *shear_damage, struct part *restrict p, const int mat_id, const float density, const float u) {

  /* Return shear_damage = 1 if already at max shear damage. */
  const float shear_damage_prev = *shear_damage;
  if (shear_damage_prev >= 1.f) {
    *shear_damage = 1.f;
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

  /* Calculate pressure and set invariant of total plastric strain. */
  const float pressure =
        gas_pressure_from_internal_energy(density, u, mat_id);
  const float strain_rate_invariant = p->strength_data.total_plastic_strain;

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
  const float shear_damage_new = fminf(strain_rate_invariant / plastic_strain_at_failure, 1.f);

  // ### Main questions here are:
  // ### 1) what the eqn is in the "else" above:
  // ### 2) Can damage decrease if e.g.  plastic_strain_at_failure increases but evolved plastic strain is const

  /* Prevent damage from decreasing. */
  *shear_damage = fmaxf(shear_damage_new, shear_damage_prev);
}

#endif /* SWIFT_DAMAGE_SHEAR_COLLINS04_H */
