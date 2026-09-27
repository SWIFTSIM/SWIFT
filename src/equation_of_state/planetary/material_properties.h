/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016   Matthieu Schaller (schaller@strw.leidenuniv.nl).
 *               2025   Jacob Kegerreis (j.kegerreis@imperial.ac.uk).
 *               2025   Thomas Sandnes (thomas.d.sandnes@durham.ac.uk).
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
#ifndef SWIFT_PLANETARY_MATERIAL_PROPS_H
#define SWIFT_PLANETARY_MATERIAL_PROPS_H

/**
 * @file equation_of_state/planetary/material_properties.h
 *
 * Material properties beyond the base EoS.
 */

/* Some standard headers. */
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "common_io.h"
#include "eos_setup.h"
#include "eos_utilities.h"
#include "error.h"
#include "inline.h"
#include "parser.h"
#include "physical_constants.h"
#include "restart.h"
#include "units.h"

/**
 * @brief Returns whether or not the material is in a solid phase and not fluid.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
material_phase_from_internal_energy(
    const float density, const float u, const enum eos_planetary_material_id mat_id) {

  const enum eos_planetary_type_id type =
      (enum eos_planetary_type_id)(mat_id / eos_type_factor);
  const int unit_id = mat_id % eos_type_factor;

  const int mat_index = material_index_from_mat_id(mat_id);

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_type_idg:
      return idg_phase_from_internal_energy(
          density, u, &eos.all_mat_params[mat_index], &eos.all_idg[unit_id]);

    /* Tillotson EoS */
    case eos_type_Til:
      return Til_phase_from_internal_energy(
          density, u, &eos.all_mat_params[mat_index], &eos.all_Til[unit_id]);

    /* Custom user-provided Tillotson EoS */
    case eos_type_Til_custom:
      return Til_phase_from_internal_energy(
          density, u, &eos.all_mat_params[mat_index],
          &eos.all_Til_custom[unit_id]);

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_type_HM80:
      return HM80_phase_from_internal_energy(
          density, u, &eos.all_mat_params[mat_index], &eos.all_HM80[unit_id]);

    /* SESAME EoS */
    case eos_type_SESAME:
      return SESAME_phase_from_internal_energy(
          density, u, &eos.all_mat_params[mat_index], &eos.all_SESAME[unit_id]);

    /* ANEOS -- using SESAME-style tables */
    case eos_type_ANEOS:
      return SESAME_phase_from_internal_energy(
          density, u, &eos.all_mat_params[mat_index], &eos.all_ANEOS[unit_id]);

    /*! Linear EoS -- user-provided parameters */
    case eos_type_linear:
      return linear_phase_from_internal_energy(
          density, u, &eos.all_mat_params[mat_index], &eos.all_linear[unit_id]);

    /*! Generic user-provided custom tables */
    case eos_type_custom:
      return SESAME_phase_from_internal_energy(
          density, u, &eos.all_mat_params[mat_index], &eos.all_custom[unit_id]);

    default:
      return -1.f;
  }
}

#ifdef MATERIAL_STRENGTH

// Material parameters

/** @brief Returns the shear modulus of a material */
__attribute__((always_inline)) INLINE static float material_shear_mod(
    const enum eos_planetary_material_id mat_id) {
  const int mat_index = material_index_from_mat_id(mat_id);
  return eos.all_mat_params[mat_index].shear_mod;
}

/** @brief Returns the bulk modulus of a material */
__attribute__((always_inline)) INLINE static float material_bulk_mod(
    const enum eos_planetary_material_id mat_id) {
  const int mat_index = material_index_from_mat_id(mat_id);
  return eos.all_mat_params[mat_index].bulk_mod;
}

/** @brief Returns the melting temperature of a material */
__attribute__((always_inline)) INLINE static float material_T_melt(
    const enum eos_planetary_material_id mat_id) {
  const int mat_index = material_index_from_mat_id(mat_id);
  return eos.all_mat_params[mat_index].T_melt;
}

/** @brief Returns the rho_0 of a material */
__attribute__((always_inline)) INLINE static float material_rho_0(
    const enum eos_planetary_material_id mat_id) {
  const int mat_index = material_index_from_mat_id(mat_id);
  return eos.all_mat_params[mat_index].rho_0;
}

#if defined(STRENGTH_YIELD_STRESS_BENZ_ASPHAUG)
  /** @brief Returns the Y_0 of a material */
  __attribute__((always_inline)) INLINE static float material_Y_0(
      const enum eos_planetary_material_id mat_id) {
    const int mat_index = material_index_from_mat_id(mat_id);
    return eos.all_mat_params[mat_index].Y_0;
  }
#elif defined(STRENGTH_YIELD_STRESS_COLLINS)
  /** @brief Returns the Y_0 of a material */
  __attribute__((always_inline)) INLINE static float material_Y_0(
      const enum eos_planetary_material_id mat_id) {
    const int mat_index = material_index_from_mat_id(mat_id);
    return eos.all_mat_params[mat_index].Y_0;
  }

  /** @brief Returns the Y_M of a material */
  __attribute__((always_inline)) INLINE static float material_Y_M(
      const enum eos_planetary_material_id mat_id) {
    const int mat_index = material_index_from_mat_id(mat_id);
    return eos.all_mat_params[mat_index].Y_M;
  }

  /** @brief Returns the mu_i of a material */
  __attribute__((always_inline)) INLINE static float material_mu_i(
      const enum eos_planetary_material_id mat_id) {
    const int mat_index = material_index_from_mat_id(mat_id);
    return eos.all_mat_params[mat_index].mu_i;
  }

  /** @brief Returns the mu_d of a material */
  __attribute__((always_inline)) INLINE static float material_mu_d(
      const enum eos_planetary_material_id mat_id) {
    const int mat_index = material_index_from_mat_id(mat_id);
    return eos.all_mat_params[mat_index].mu_d;
  }
#endif

#if defined(STRENGTH_DAMAGE_SHEAR_COLLINS)
  /** @brief Returns the brittle to ductile transition pressure of a material */
  __attribute__((always_inline)) INLINE static float
  material_brittle_to_ductile_pressure(const enum eos_planetary_material_id mat_id) {
    const int mat_index = material_index_from_mat_id(mat_id);
    return eos.all_mat_params[mat_index].brittle_to_ductile_pressure;
  }

  /** @brief Returns the brittle to plastic transition pressure of a material */
  __attribute__((always_inline)) INLINE static float
  material_brittle_to_plastic_pressure(const enum eos_planetary_material_id mat_id) {
    const int mat_index = material_index_from_mat_id(mat_id);
    return eos.all_mat_params[mat_index].brittle_to_plastic_pressure;
  }
#endif

// Method-specific material parameters

#if defined(STRENGTH_YIELD_STRESS_WEAKENING_THERMAL)
  /** @brief Returns the yield stress thermal weakening parameter of a material */
  __attribute__((always_inline)) INLINE static float
  material_yield_weakening_thermal_xi(const enum eos_planetary_material_id mat_id) {
    const int mat_index = material_index_from_mat_id(mat_id);
    return eos.all_mat_params[mat_index].yield_weakening_thermal_xi;
  }
#endif /* STRENGTH_YIELD_STRESS_WEAKENING_THERMAL */

#if defined(STRENGTH_YIELD_STRESS_WEAKENING_DENSITY)
  /** @brief Returns the yield stress density weakening multiplication parameter
   * of a material */
  __attribute__((always_inline)) INLINE static float
  material_yield_weakening_density_mult_param(const enum eos_planetary_material_id mat_id) {
    const int mat_index = material_index_from_mat_id(mat_id);
    return eos.all_mat_params[mat_index].yield_weakening_density_mult_param;
  }

  /** @brief Returns the yield stress density weakening exponent parameter of a
   * material */
  __attribute__((always_inline)) INLINE static float
  material_yield_weakening_density_pow_param(const enum eos_planetary_material_id mat_id) {
    const int mat_index = material_index_from_mat_id(mat_id);
    return eos.all_mat_params[mat_index].yield_weakening_density_pow_param;
  }
#endif /* STRENGTH_YIELD_STRESS_WEAKENING_DENSITY */

#if defined(STRENGTH_ARTIFICIAL_STRESS_MON2000)
  /** @brief Returns the artificial stress n parameter of a material */
  __attribute__((always_inline)) INLINE static float material_artif_stress_n(
      const enum eos_planetary_material_id mat_id) {
    const int mat_index = material_index_from_mat_id(mat_id);
    return eos.all_mat_params[mat_index].artif_stress_n;
  }

  /** @brief Returns the artificial stress epsilon parameter of a material */
  __attribute__((always_inline)) INLINE static float material_artif_stress_epsilon(
      const enum eos_planetary_material_id mat_id) {
    const int mat_index = material_index_from_mat_id(mat_id);
    return eos.all_mat_params[mat_index].artif_stress_epsilon;
  }
#endif /* STRENGTH_ARTIFICIAL_STRESS_MON2000 */

#endif /* MATERIAL_STRENGTH */

/**
 * @brief Set the default material parameters (SI units).
 *
 * Any of these are overwritten by parameters given in a material parameter file.
 *
 * @param mat_params The material parameters.
 */
INLINE static void set_material_params_default(struct mat_params *mat_params) {

  // General properties
  mat_params->state_type = mat_state_type_fluid;

#ifdef MATERIAL_STRENGTH
  // General material strength parameters
  mat_params->shear_mod = 0.f;
  mat_params->bulk_mod = 0.f;
  mat_params->T_melt = 0.f;
  mat_params->rho_0 = 0.f;

  // Specific constants for material-strength schemes
  #if defined(STRENGTH_YIELD_STRESS_BENZ_ASPHAUG)
    mat_params->Y_0 = 0.f;
  #elif defined(STRENGTH_YIELD_STRESS_COLLINS)
    mat_params->Y_0 = 0.f;
    mat_params->Y_M = 0.f;
    mat_params->mu_i = 0.f;
    mat_params->mu_d = 0.f;
  #endif

  #if defined(STRENGTH_YIELD_STRESS_WEAKENING_THERMAL)
    mat_params->yield_weakening_thermal_xi = 1.2f;
  #endif

  #if defined(STRENGTH_YIELD_STRESS_WEAKENING_DENSITY)
    mat_params->yield_weakening_density_mult_param = 0.85f;
    mat_params->yield_weakening_density_pow_param = 4.f;
  #endif

  #if defined(STRENGTH_DAMAGE_SHEAR_COLLINS)
    mat_params->brittle_to_ductile_pressure = 0.f;
    mat_params->brittle_to_plastic_pressure = 0.f;
  #endif

  #if defined(STRENGTH_ARTIFICIAL_STRESS_MON2000)
    mat_params->artif_stress_n = 4.f;
    mat_params->artif_stress_epsilon = 0.2f;
  #endif
#endif /* MATERIAL_STRENGTH */
}

/**
 * @brief Overwrite material parameters with any given in a material parameter file.
 *
 * Parameters that are not in the file keep their current values.
 *
 * @param mat_params The material parameters.
 * @param param_file The material parameter file.
 */
INLINE static void set_material_params_from_file(struct mat_params *mat_params,
                                                 char *param_file) {

#ifdef MATERIAL_STRENGTH
  // Load parameter file
  struct swift_params *file_params =
      (struct swift_params *)malloc(sizeof(struct swift_params));
  parser_read_file(param_file, file_params);

  // General properties
  mat_params->state_type = (enum mat_state_type)parser_get_opt_param_int(
      file_params, "Material:state_type", mat_params->state_type);

  // General material strength parameters
  mat_params->shear_mod = parser_get_opt_param_float(
      file_params, "Strength:shear_mod", mat_params->shear_mod);
  mat_params->bulk_mod = parser_get_opt_param_float(
      file_params, "Strength:bulk_mod", mat_params->bulk_mod);
  mat_params->T_melt = parser_get_opt_param_float(
      file_params, "Strength:T_melt", mat_params->T_melt);
  mat_params->rho_0 = parser_get_opt_param_float(
      file_params, "Strength:rho_0", mat_params->rho_0);

  // Specific constants for material-strength schemes
  #if defined(STRENGTH_YIELD_STRESS_BENZ_ASPHAUG)
    mat_params->Y_0 = parser_get_opt_param_float(
        file_params, "YieldStress:Y_0", mat_params->Y_0);
  #elif defined(STRENGTH_YIELD_STRESS_COLLINS)
    mat_params->Y_0 = parser_get_opt_param_float(
        file_params, "YieldStress:Y_0", mat_params->Y_0);
    mat_params->Y_M = parser_get_opt_param_float(
        file_params, "YieldStress:Y_M", mat_params->Y_M);
    mat_params->mu_i = parser_get_opt_param_float(
        file_params, "YieldStress:mu_i", mat_params->mu_i);
    mat_params->mu_d = parser_get_opt_param_float(
        file_params, "YieldStress:mu_d", mat_params->mu_d);
  #endif

  #if defined(STRENGTH_YIELD_STRESS_WEAKENING_THERMAL)
    mat_params->yield_weakening_thermal_xi = parser_get_opt_param_float(
        file_params, "YieldStress:yield_weakening_thermal_xi",
        mat_params->yield_weakening_thermal_xi);
  #endif

  #if defined(STRENGTH_YIELD_STRESS_WEAKENING_DENSITY)
    mat_params->yield_weakening_density_mult_param = parser_get_opt_param_float(
        file_params, "YieldStress:yield_weakening_density_mult_param",
        mat_params->yield_weakening_density_mult_param);
    mat_params->yield_weakening_density_pow_param = parser_get_opt_param_float(
        file_params, "YieldStress:yield_weakening_density_pow_param",
        mat_params->yield_weakening_density_pow_param);
  #endif

  #if defined(STRENGTH_DAMAGE_SHEAR_COLLINS)
    mat_params->brittle_to_ductile_pressure = parser_get_opt_param_float(
        file_params, "DamageShearCollins:brittle_to_ductile_pressure",
        mat_params->brittle_to_ductile_pressure);
    mat_params->brittle_to_plastic_pressure = parser_get_opt_param_float(
        file_params, "DamageShearCollins:brittle_to_plastic_pressure",
        mat_params->brittle_to_plastic_pressure);
  #endif

  #if defined(STRENGTH_ARTIFICIAL_STRESS_MON2000)
    mat_params->artif_stress_n = parser_get_opt_param_float(
        file_params, "ArtificialStress:n", mat_params->artif_stress_n);
    mat_params->artif_stress_epsilon = parser_get_opt_param_float(
        file_params, "ArtificialStress:epsilon", mat_params->artif_stress_epsilon);
  #endif

  free(file_params);
#endif /* MATERIAL_STRENGTH */
}

/**
 * @brief Convert units of material parameters from SI to internal units.
 *
 * @param mat_params The material parameters.
 * @param us The internal unit system.
 */
INLINE static void convert_units_material_params(struct mat_params *mat_params,
                                                 const struct unit_system *us) {

#ifdef MATERIAL_STRENGTH
  struct unit_system si;
  units_init_si(&si);

  // General material strength parameters
  // SI to cgs
  mat_params->shear_mod *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);
  mat_params->bulk_mod *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);
  mat_params->T_melt *= units_cgs_conversion_factor(&si, UNIT_CONV_TEMPERATURE);
  mat_params->rho_0 *= units_cgs_conversion_factor(&si, UNIT_CONV_DENSITY);

  // cgs to internal
  mat_params->shear_mod /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
  mat_params->bulk_mod /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
  mat_params->T_melt /= units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
  mat_params->rho_0 /= units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);

  // Specific constants for material-strength schemes
  #if defined(STRENGTH_YIELD_STRESS_BENZ_ASPHAUG)
    // SI to cgs
    mat_params->Y_0 *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);

    // cgs to internal
    mat_params->Y_0 /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
  #elif defined(STRENGTH_YIELD_STRESS_COLLINS)
    // SI to cgs
    mat_params->Y_0 *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);
    mat_params->Y_M *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);

    // cgs to internal
    mat_params->Y_0 /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
    mat_params->Y_M /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
  #endif

  #if defined(STRENGTH_DAMAGE_SHEAR_COLLINS)
    // SI to cgs
    mat_params->brittle_to_ductile_pressure *=
        units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);
    mat_params->brittle_to_plastic_pressure *=
        units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);

    // cgs to internal
    mat_params->brittle_to_ductile_pressure /=
        units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
    mat_params->brittle_to_plastic_pressure /=
        units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
  #endif
#endif /* MATERIAL_STRENGTH */
}

/**
 * @brief Check material parameters that would otherwise silently remove
 * strength or give divisions by zero, for materials that can be solid.
 *
 * @param mat_params The material parameters.
 * @param mat_id The material ID.
 */
INLINE static void check_material_params(const struct mat_params *mat_params,
                                         const enum eos_planetary_material_id mat_id) {

#ifdef MATERIAL_STRENGTH
  if (mat_params->state_type == mat_state_type_fluid) return;

  #if defined(STRENGTH_YIELD_STRESS_WEAKENING_THERMAL)
    if ((mat_params->T_melt <= 0.f) ||
        (mat_params->yield_weakening_thermal_xi <= 0.f))
      error("Material %d: thermal yield weakening needs Strength:T_melt > 0 "
            "and YieldStress:yield_weakening_thermal_xi > 0", mat_id);
  #endif

  #if defined(STRENGTH_YIELD_STRESS_WEAKENING_DENSITY)
    if (mat_params->rho_0 <= 0.f)
      error("Material %d: density yield weakening needs Strength:rho_0 > 0",
            mat_id);
  #endif

  #if defined(STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG)
    if ((mat_params->bulk_mod <= 0.f) || (mat_params->shear_mod <= 0.f))
      error("Material %d: tensile damage needs Strength:bulk_mod > 0 and "
            "Strength:shear_mod > 0", mat_id);
  #endif

  #if defined(STRENGTH_DAMAGE_SHEAR_COLLINS)
    if ((mat_params->brittle_to_ductile_pressure <= 0.f) ||
        (mat_params->brittle_to_plastic_pressure <=
         mat_params->brittle_to_ductile_pressure))
      error("Material %d: shear damage needs 0 < "
            "DamageShearCollins:brittle_to_ductile_pressure < "
            "DamageShearCollins:brittle_to_plastic_pressure", mat_id);
  #endif

  #if defined(STRENGTH_ARTIFICIAL_STRESS_MON2000)
    if ((mat_params->artif_stress_epsilon > 0.f) &&
        (mat_params->artif_stress_n <= 0.f))
      error("Material %d: ArtificialStress:n must be > 0", mat_id);
  #endif
#endif /* MATERIAL_STRENGTH */
}

/**
 * @brief Set the material parameters of a material.
 *
 * Starts from the default parameters, which are overwritten by any parameters
 * given in the material parameter file, if one is provided.
 *
 * @param all_mat_params The material parameters of all materials.
 * @param mat_id The material ID.
 * @param param_file The material parameter file, or "NoFile".
 * @param us The internal unit system.
 */
INLINE static void set_material_params(struct mat_params *all_mat_params,
                                       enum eos_planetary_material_id mat_id,
                                       char *param_file,
                                       const struct unit_system *us) {

  const int mat_index = material_index_from_mat_id(mat_id);
  struct mat_params *mat_params = &all_mat_params[mat_index];

  // Default parameters
  set_material_params_default(mat_params);

  // Overwrite with any parameters given in the material parameter file
  if (strcmp(param_file, "NoFile") != 0) {
    set_material_params_from_file(mat_params, param_file);
  }

  // Convert units
  convert_units_material_params(mat_params, us);

  // Check parameters
  check_material_params(mat_params, mat_id);
}

#endif /* SWIFT_PLANETARY_MATERIAL_PROPS_H */
