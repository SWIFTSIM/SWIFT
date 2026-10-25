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
#ifndef SWIFT_WEAKENING_DENSITY_H
#define SWIFT_WEAKENING_DENSITY_H

/**
 * @file strength/yield_stress_weakening/weakening_density.h
 * @brief Density weakening of the yield stress.
 */

#include "math.h"

/**
 * @brief Sets the density weakening properties of the strength state from the particle.
 *
 * Nothing extra is needed in the state for this method.
 *
 * @param s The strength state of the particle.
 * @param p The particle of interest.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void yield_weakening_density_set_state(
    struct strength_state *s, const struct part *restrict p,
    const struct mat_params *mp) {}

/**
 * @brief Apply density weakening to the yield stress.
 *
 * Weaken the material's yield stress as its density decreases below a certain
 * point. This method is discussed by Luther+2022, where it is described in the
 * context of its implementation in iSALE.
 *
 * Method parameters needed in material parameter file:
 * Strength:
 *     rho_0: Material rho_0.
 * YieldWeakeningDensity:
 *     mult_param: Density weakening muliplication factor (default 0.85).
 *     pow_param: Density weakening power factor (default 4).
 *
 * @param Y The yield stress to be weakened.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void
yield_weakening_apply_density_to_yield_stress(float *Y,
                                              struct strength_state *s,
                                              const struct mat_params *mp) {

  /* Method parameters. */
  const float rho_0 = mp->rho_0;
  const float a = mp->yield_weakening_density_mult_param;
  const float rho_weak = a * rho_0;
  const float density = s->density;

  /* Only apply weakening below a density threshold. */
  if (density < rho_weak) {
    const float b = mp->yield_weakening_density_pow_param;

    /* Apply weakening. */
    *Y *= powf(density / rho_weak, b);
  }
}

#endif /* SWIFT_WEAKENING_DENSITY_H */
