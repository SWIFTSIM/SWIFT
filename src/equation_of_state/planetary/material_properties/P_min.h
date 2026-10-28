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
 * pressure at which the intact yield stress reaches zero in tension.
 *
 * Only for the Collins yield stress with mu_i > 0 and Y_M > Y_0, since
 * otherwise the yield stress does not depend on pressure.
 *
 * @param mp The material parameters.
 * @param mat_id The material ID.
 */
__attribute__((always_inline)) INLINE static float material_P_min_yield_envelope(
    const struct mat_params *mp, const enum eos_planetary_material_id mat_id) {

#if defined(STRENGTH_YIELD_STRESS_COLLINS)
  const float Y_0 = mp->yield_stress_collins_Y_0;
  const float Y_M = mp->yield_stress_collins_Y_M;
  const float mu_i = mp->yield_stress_collins_mu_i;
  if ((mu_i <= 0.f) || (Y_M <= Y_0))
    error("Material %d: P_min_scheme yield_envelope needs "
          "YieldStressCollins:mu_i > 0 and Y_M > Y_0", mat_id);

  return -Y_0 * (Y_M - Y_0) / (mu_i * Y_M);
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
