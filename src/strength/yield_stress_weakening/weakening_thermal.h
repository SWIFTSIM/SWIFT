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
#ifndef SWIFT_WEAKENING_THERMAL_H
#define SWIFT_WEAKENING_THERMAL_H

/**
 * @file strength/yield_stress_weakening/weakening_thermal.h
 * @brief Thermal weakening of the yield stress.
 */

#include "math.h"

/**
 * @brief Sets the thermal weakening properties of the strength state from the particle.
 *
 * Nothing extra is needed in the state for this method.
 *
 * @param s The strength state of the particle.
 * @param p The particle of interest.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void yield_weakening_thermal_set_state(
    struct strength_state *s, const struct part *restrict p,
    const struct mat_params *mp) {}

/**
 * @brief Apply thermal weakening to the yield stress.
 *
 * Weaken the material's yield stress as its temperature approaches the
 * melting temperature. This method is presented by Ohnaka1995, and
 * is used by e.g. Collins+2004 and Emsenhuber+2018.
 *
 * Method parameters needed in material parameter file:
 * Strength:
 *     T_melt: Melting temperature (K).
 * YieldStress:
 *     yield_weakening_thermal_xi: Thermal weakening xi factor (default 1.2).
 *
 * @param Y The yield stress to be weakened.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void
yield_weakening_apply_temperature_to_yield_stress(float *Y,
                                                  const struct strength_state *s,
                                                  const struct mat_params *mp) {

  /* Method parameters. */
  const float xi = mp->yield_weakening_thermal_xi;
  const float T_melt = mp->T_melt;
  const float temperature = s->temperature;

  /* Apply weakening. */
  if (temperature > 0.f && temperature < T_melt) {
    *Y *= tanhf(xi * (T_melt / temperature - 1.f));
  } else if (temperature >= T_melt) {
    *Y = 0.f;
  }
}

#endif /* SWIFT_WEAKENING_THERMAL_H */
