/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
 *               2024 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
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
#ifndef SWIFT_PLANETARY_EOS_UTILITIES_H
#define SWIFT_PLANETARY_EOS_UTILITIES_H

/**
 * @file equation_of_state/planetary/eos_utilities.h
 *
 * Utitilies for the planetary equations of state.
 */

/* Local headers. */
#include "eos_hm80.h"
#include "eos_ideal_gas.h"
#include "eos_linear.h"
#include "eos_sesame.h"
#include "eos_setup.h"
#include "eos_tillotson.h"
#include "parser.h"

/**
 * @brief The parameters of the equation of state.
 *
 * Parameter structs for all EoS of each type.
 */
struct eos_parameters {
  // Primary EoS parameters, with different elements for each EoS type
  struct idg_params all_idg[eos_count_idg];
  struct Til_params all_Til[eos_count_Til];
  struct Til_params all_Til_custom[eos_count_Til_custom];
  struct HM80_params all_HM80[eos_count_HM80];
  struct SESAME_params all_SESAME[eos_count_SESAME];
  struct SESAME_params all_ANEOS[eos_count_ANEOS];
  struct linear_params all_linear[eos_count_linear];
  struct SESAME_params all_custom[eos_count_custom];

  // Material parameters, for any EoS type
  struct mat_params all_mat_params[eos_count_total];
  // Method parameters, same for all EoS type
  struct method_params method_params;
};

/*! Primary EoS parameter struct */
extern struct eos_parameters eos;

/**
 * @brief Returns whether or not the material is in a solid state and not fluid.
 */
__attribute__((always_inline)) INLINE static int material_index_from_mat_id(
    enum eos_planetary_material_id mat_id) {

  const enum eos_planetary_type_id type =
      (enum eos_planetary_type_id)(mat_id / eos_type_factor);
  const int unit_id = mat_id % eos_type_factor;

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_type_idg:
      return eos_cumul_count_idg + unit_id;

    /* Tillotson EoS */
    case eos_type_Til:
      return eos_cumul_count_Til + unit_id;

    /* Custom user-provided Tillotson EoS */
    case eos_type_Til_custom:
      return eos_cumul_count_Til_custom + unit_id;

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_type_HM80:
      return eos_cumul_count_HM80 + unit_id;

    /* SESAME EoS */
    case eos_type_SESAME:
      return eos_cumul_count_SESAME + unit_id;

    /* ANEOS -- using SESAME-style tables */
    case eos_type_ANEOS:
      return eos_cumul_count_ANEOS + unit_id;

    /*! Linear EoS -- user-provided parameters */
    case eos_type_linear:
      return eos_cumul_count_linear + unit_id;

    /*! Generic user-provided custom tables */
    case eos_type_custom:
      return eos_cumul_count_custom + unit_id;

    default:
      return -1.f;
  }
}

/**
 * @brief Convert units of material parameters from SI to internal units.
 */
INLINE static void convert_units_material_params(struct mat_params *mat_params,
                                       struct method_params *method_params,
                                       const struct unit_system *us) {

   struct unit_system si;
   units_init_si(&si);

#ifdef MATERIAL_STRENGTH
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

   // Specific EoS-dependent constants for material-strength schemes
   #if defined(STRENGTH_YIELD_BENZ_ASPHAUG)
     mat_params->Y_0 *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);

     // cgs to internal
     mat_params->Y_0 /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
   #elif defined(STRENGTH_YIELD_COLLINS)
     // SI to cgs
     mat_params->Y_0 *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);
     mat_params->Y_M *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);

     // cgs to internal
     mat_params->Y_0 /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
     mat_params->Y_M /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
   #endif

   #if defined(STRENGTH_DAMAGE_SHEAR_COLLINS)
     // SI to cgs
     mat_params->brittle_to_ductile_pressure *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);
     mat_params->brittle_to_plastic_pressure *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);

     // cgs to internal
     mat_params->brittle_to_ductile_pressure /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
     mat_params->brittle_to_plastic_pressure /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
   #endif
#endif /* MATERIAL_STRENGTH */

 }
/**
 * @brief Load material parameters from a file.
 */
INLINE static void set_material_params(struct mat_params *all_mat_params,
                                       struct method_params *method_params,
                                       enum eos_planetary_material_id mat_id,
                                       char *param_file,
                                       const struct unit_system *us) {

  const int mat_index = material_index_from_mat_id(mat_id);
  struct mat_params *mat_params = &all_mat_params[mat_index];

#ifdef MATERIAL_STRENGTH
  // Load parameter file
  struct swift_params *file_params =
      (struct swift_params *)malloc(sizeof(struct swift_params));
  parser_read_file(param_file, file_params);

  // General properties
  int phase_state = parser_get_opt_param_int(
      file_params, "Material:phase_state", mat_phase_state_fluid);
  mat_params->phase_state = (enum mat_phase_state)phase_state;

  // General material strength parameters
  mat_params->shear_mod =
      parser_get_opt_param_float(file_params, "Strength:shear_mod", 0.f);
  mat_params->bulk_mod =
      parser_get_opt_param_float(file_params, "Strength:bulk_mod", 0.f);
  mat_params->T_melt =
      parser_get_opt_param_float(file_params, "Strength:T_melt", 0.f);
  mat_params->rho_0 =
      parser_get_opt_param_float(file_params, "Strength:rho_0", 0.f);


  // Specific EoS-independent constants for material-strength schemes
  #if defined(STRENGTH_STRESS_MON2000) || defined(STRENGTH_STRESS_BASIS_INDP)
    method_params->artif_stress_n = parser_get_opt_param_float(
        file_params, "ArtificialStress:n", 0.f);
    method_params->artif_stress_epsilon = parser_get_opt_param_float(
        file_params, "ArtificialStress:epsilon", 0.f);
  #endif /* STRENGTH_STRESS_MON2000 || STRENGTH_STRESS_BASIS_INDP */

  // Specific EoS-dependent constants for material-strength schemes
  #if defined(STRENGTH_YIELD_BENZ_ASPHAUG)
    mat_params->Y_0 =
        parser_get_opt_param_float(file_params, "YieldStress:Y_0", 0.f);
  #elif defined(STRENGTH_YIELD_COLLINS)
    mat_params->Y_0 =
        parser_get_opt_param_float(file_params, "YieldStress:Y_0", 0.f);
    mat_params->Y_M =
        parser_get_opt_param_float(file_params, "YieldStress:Y_M", 0.f);
    mat_params->mu_i =
        parser_get_opt_param_float(file_params, "YieldStress:mu_i", 0.f);
    mat_params->mu_d =
        parser_get_opt_param_float(file_params, "YieldStress:mu_d", 0.f);
  #endif

  #if defined(STRENGTH_YIELD_THERMAL_SOFTENING)
    method_params->yield_thermal_soft_xi = parser_get_opt_param_float(
        file_params, "YieldStress:yield_thermal_soft_xi", 0.f);
  #endif

  #if defined(STRENGTH_YIELD_DENSITY_SOFTENING)
    method_params->yield_density_soft_mult_param = parser_get_opt_param_float(
        file_params, "YieldStress:yield_density_soft_mult_param", 0.f);
    method_params->yield_density_soft_pow_param = parser_get_opt_param_float(
        file_params, "YieldStress:yield_density_soft_pow_param", 0.f);
  #endif

  #if defined(STRENGTH_DAMAGE_SHEAR_COLLINS)
    mat_params->brittle_to_ductile_pressure = parser_get_opt_param_float(
        file_params, "DamageShearCollins:brittle_to_ductile_pressure", 0.f);
    mat_params->brittle_to_plastic_pressure = parser_get_opt_param_float(
        file_params, "DamageShearCollins:brittle_to_plastic_pressure", 0.f);
  #endif

#else
  mat_params->phase_state = mat_phase_state_fluid;
#endif /* MATERIAL_STRENGTH */

  // Convert units
  convert_units_material_params(mat_params, method_params, us);

}

#endif /* SWIFT_PLANETARY_EOS_UTILITIES_H */
