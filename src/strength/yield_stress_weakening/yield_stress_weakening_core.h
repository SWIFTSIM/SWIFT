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
#ifndef SWIFT_YIELD_STRESS_WEAKENING_CORE_H
#define SWIFT_YIELD_STRESS_WEAKENING_CORE_H

/**
 * @file strength/yield_stress_weakening/yield_stress_weakening_core.h
 * @brief Core yield stress weakening scheme.
 */

/**
 * @brief Sets the yield stress weakening properties of the strength state from
 * the particle.
 *
 * @param s The strength state of the particle.
 * @param p The particle of interest.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void yield_weakening_set_state(
    struct strength_state *s, const struct part *restrict p,
    const struct mat_params *mp) {

  yield_weakening_thermal_set_state(s, p, mp);
  yield_weakening_density_set_state(s, p, mp);
}

/**
 * @brief Apply the thermal and density weakening to the yield stress.
 *
 * @param Y The yield stress to be weakened.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void
yield_weakening_apply_to_yield_stress(float *Y, struct strength_state *s,
                                      const struct mat_params *mp) {

  yield_weakening_apply_density_to_yield_stress(Y, s, mp);
  yield_weakening_apply_temperature_to_yield_stress(Y, s, mp);
}

#endif /* SWIFT_YIELD_STRESS_WEAKENING_CORE_H */
