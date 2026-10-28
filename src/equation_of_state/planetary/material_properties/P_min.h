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
#ifndef SWIFT_PLANETARY_MATERIAL_PROPERTIES_P_MIN_H
#define SWIFT_PLANETARY_MATERIAL_PROPERTIES_P_MIN_H

/**
 * @file equation_of_state/planetary/material_properties/P_min.h
 * @brief Schemes for the minimum pressure of a material.
 */

#include <float.h>
#include <math.h>

#include "../eos_setup.h"
#include "error.h"
#include "parser.h"
#include "units.h"

#ifdef MATERIAL_STRENGTH

/**
 * @brief The constant scheme for the minimum pressure.
 *
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float material_P_min_constant(
    const struct mat_params *mp) {

  return mp->P_min_constant;
}

/**
 * @brief The yield-envelope scheme for the minimum pressure, which uses the
 * pressure at which the yield stress reaches zero in tension.
 *
 * @param mp The material parameters.
 * @param damage The total damage.
 * @param mat_id The material ID.
 */
__attribute__((always_inline)) INLINE static float material_P_min_yield_envelope(
    const struct mat_params *mp, const float damage,
    const enum eos_planetary_material_id mat_id) {

#if defined(STRENGTH_YIELD_STRESS_COLLINS)
  const float Y_0i = mp->yield_stress_collins_Y_0i;
  const float Y_Mi = mp->yield_stress_collins_Y_Mi;
  const float mu_i = mp->yield_stress_collins_mu_i;
  const float Y_0d = mp->yield_stress_collins_Y_0d;
  const float mu_d = mp->yield_stress_collins_mu_d;
  const float D = damage;

  if (mu_d <= 0.f)
    error("Material %d: P_min_scheme yield_envelope needs "
          "YieldStressCollins:mu_d > 0", mat_id);

  /* The yield stress is Y = (1 - D) Y_i + D Y_d, with Y_d not capped at zero.
   * Find the pressure at which Y is zero. */

  /* Constant Y_i = Y_0i, when it does not depend on pressure. Intact material
   * then never reaches zero, and otherwise the zero is linear in P. */
  if ((Y_Mi <= Y_0i) || (mu_i == 0.f)) {
    if (D <= 0.f) {
      return -FLT_MAX;
    }
    return -((1.f - D) * Y_0i + D * Y_0d) / (D * mu_d);
  }

  /* Solve quadratic equation for P. */
  const float k = mu_i / (Y_Mi - Y_0i);
  const float a = D * mu_d * k;
  const float b = (1.f - D) * (mu_i + Y_0i * k) + D * (mu_d + Y_0d * k);
  const float c = (1.f - D) * Y_0i + D * Y_0d;

  /* Use the larger root, the first zero reached in tension down from P = 0:
   * (-b + sqrt(b^2 - 4 a c)) / (2 a), rewritten to avoid dividing by a. */
  return -2.f * c / (b + sqrtf(fmaxf(b * b - 4.f * a * c, 0.f)));
#else
  error("Material %d: P_min_scheme yield_envelope needs the Collins yield "
        "stress", mat_id);
  return 0.f;
#endif
}

/**
 * @brief Set the default minimum pressure parameters.
 *
 * @param mp The material parameters.
 */
INLINE static void material_P_min_set_default(struct mat_params *mp) {

  mp->P_min_scheme = P_min_scheme_none;
  mp->P_min_constant = 0.f;
}

/**
 * @brief Read the minimum pressure parameters from a material parameter file.
 *
 * Only the parameters of the chosen scheme are read, and they are compulsory.
 *
 * @param mp The material parameters.
 * @param file_params The material parameter file.
 * @param param_file The name of the material parameter file.
 */
INLINE static void material_P_min_read(struct mat_params *mp,
                                       struct swift_params *file_params,
                                       const char *param_file) {

  // Scheme
  mp->P_min_scheme = (enum P_min_scheme)parser_get_opt_param_int(
      file_params, "Material:P_min_scheme", mp->P_min_scheme);

  // Specific constants for the scheme
  switch (mp->P_min_scheme) {
    case P_min_scheme_none:
      break;

    case P_min_scheme_constant:
      mp->P_min_constant = parser_get_param_float(file_params, "Material:P_min_constant");
      break;

    case P_min_scheme_yield_envelope:
      break;

    default:
      error("%s: unknown Material:P_min_scheme %d", param_file, mp->P_min_scheme);
  }
}

/**
 * @brief Convert the minimum pressure parameters from SI to internal units.
 *
 * @param mp The material parameters.
 * @param us The internal unit system.
 */
INLINE static void material_P_min_convert_units(struct mat_params *mp,
                                                const struct unit_system *us) {

  struct unit_system si;
  units_init_si(&si);

  mp->P_min_constant *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);
  mp->P_min_constant /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
}

#endif /* MATERIAL_STRENGTH */

#endif /* SWIFT_PLANETARY_MATERIAL_PROPERTIES_P_MIN_H */
