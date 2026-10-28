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
#ifndef SWIFT_PLANETARY_MATERIAL_PROPERTIES_RHO_0_H
#define SWIFT_PLANETARY_MATERIAL_PROPERTIES_RHO_0_H

/**
 * @file equation_of_state/planetary/material_properties/rho_0.h
 * @brief Schemes for the reference density of a material.
 */

#include "../eos_setup.h"
#include "../eos_utilities.h"
#include "error.h"
#include "parser.h"
#include "units.h"

#ifdef MATERIAL_STRENGTH

/**
 * @brief The constant scheme for the reference density.
 *
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float material_rho_0_constant(
    const struct mat_params *mp) {

  return mp->rho_0_constant;
}

/**
 * @brief The from-EoS scheme for the reference density, which uses the rho_0
 * parameter of the material's EoS.
 *
 * @param mat_id The material ID.
 */
__attribute__((always_inline)) INLINE static float material_rho_0_from_eos(
    const enum eos_planetary_material_id mat_id) {

  const enum eos_planetary_type_id type =
      (enum eos_planetary_type_id)(mat_id / eos_type_factor);
  const int unit_id = mat_id % eos_type_factor;

  switch (type) {
    case eos_type_Til:
      return eos.all_Til[unit_id].rho_0;

    case eos_type_Til_custom:
      return eos.all_Til_custom[unit_id].rho_0;

    case eos_type_linear:
      return eos.all_linear[unit_id].rho_0;

    default:
      error("Material %d: rho_0_scheme from_eos is not compatible with EoS type %d",
            mat_id, type);
      return 0.f;
  }
}

/**
 * @brief Set the default reference density parameters.
 *
 * @param mp The material parameters.
 */
INLINE static void material_rho_0_set_default(struct mat_params *mp) {

  mp->rho_0_scheme = rho_0_scheme_none;
  mp->rho_0_constant = 0.f;
}

/**
 * @brief Read the reference density parameters from a material parameter file.
 *
 * Only the parameters of the chosen scheme are read, and they are compulsory.
 *
 * @param mp The material parameters.
 * @param file_params The material parameter file.
 * @param param_file The name of the material parameter file.
 */
INLINE static void material_rho_0_read(struct mat_params *mp,
                                       struct swift_params *file_params,
                                       const char *param_file) {

  // Scheme
  mp->rho_0_scheme = (enum rho_0_scheme)parser_get_opt_param_int(
      file_params, "Material:rho_0_scheme", mp->rho_0_scheme);

  // Specific constants for the scheme
  switch (mp->rho_0_scheme) {
    case rho_0_scheme_none:
      break;

    case rho_0_scheme_constant:
      mp->rho_0_constant = parser_get_param_float(file_params, "Material:rho_0_constant");
      if (mp->rho_0_constant <= 0.f)
        error("%s: Material:rho_0_constant must be > 0", param_file);
      break;

    case rho_0_scheme_from_eos:
      break;

    default:
      error("%s: unknown Material:rho_0_scheme %d", param_file, mp->rho_0_scheme);
  }
}

/**
 * @brief Convert the reference density parameters from SI to internal units.
 *
 * @param mp The material parameters.
 * @param us The internal unit system.
 */
INLINE static void material_rho_0_convert_units(struct mat_params *mp,
                                                const struct unit_system *us) {

  struct unit_system si;
  units_init_si(&si);

  mp->rho_0_constant *= units_cgs_conversion_factor(&si, UNIT_CONV_DENSITY);
  mp->rho_0_constant /= units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
}

#endif /* MATERIAL_STRENGTH */

#endif /* SWIFT_PLANETARY_MATERIAL_PROPERTIES_RHO_0_H */
