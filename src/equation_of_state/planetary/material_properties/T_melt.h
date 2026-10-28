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
#ifndef SWIFT_PLANETARY_MATERIAL_PROPERTIES_T_MELT_H
#define SWIFT_PLANETARY_MATERIAL_PROPERTIES_T_MELT_H

/**
 * @file equation_of_state/planetary/material_properties/T_melt.h
 * @brief Schemes for the melting temperature of a material.
 */

#include "../eos_setup.h"
#include "error.h"
#include "parser.h"
#include "units.h"

#ifdef MATERIAL_STRENGTH

/**
 * @brief The constant scheme for the melting temperature.
 *
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float material_T_melt_constant(
    const struct mat_params *mp) {

  return mp->T_melt_constant;
}

/**
 * @brief Set the default melting temperature parameters.
 *
 * @param mp The material parameters.
 */
INLINE static void material_T_melt_set_default(struct mat_params *mp) {

  mp->T_melt_scheme = T_melt_scheme_none;
  mp->T_melt_constant = 0.f;
}

/**
 * @brief Read the melting temperature parameters from a material parameter file.
 *
 * Only the parameters of the chosen scheme are read, and they are compulsory.
 *
 * @param mp The material parameters.
 * @param file_params The material parameter file.
 * @param param_file The name of the material parameter file.
 */
INLINE static void material_T_melt_read(struct mat_params *mp,
                                        struct swift_params *file_params,
                                        const char *param_file) {

  // Scheme
  mp->T_melt_scheme = (enum T_melt_scheme)parser_get_opt_param_int(
      file_params, "Material:T_melt_scheme", mp->T_melt_scheme);

  // Specific constants for the scheme
  switch (mp->T_melt_scheme) {
    case T_melt_scheme_none:
      break;

    case T_melt_scheme_constant:
      mp->T_melt_constant = parser_get_param_float(file_params, "Material:T_melt_constant");
      if (mp->T_melt_constant <= 0.f)
        error("%s: Material:T_melt_constant must be > 0", param_file);
      break;

    default:
      error("%s: unknown Material:T_melt_scheme %d", param_file, mp->T_melt_scheme);
  }
}

/**
 * @brief Convert the melting temperature parameters from SI to internal units.
 *
 * @param mp The material parameters.
 * @param us The internal unit system.
 */
INLINE static void material_T_melt_convert_units(struct mat_params *mp,
                                                 const struct unit_system *us) {

  struct unit_system si;
  units_init_si(&si);

  mp->T_melt_constant *= units_cgs_conversion_factor(&si, UNIT_CONV_TEMPERATURE);
  mp->T_melt_constant /= units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
}

#endif /* MATERIAL_STRENGTH */

#endif /* SWIFT_PLANETARY_MATERIAL_PROPERTIES_T_MELT_H */
