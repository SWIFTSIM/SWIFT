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
#include <float.h>
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

/* Schemes for material properties */
#include "material_properties/T_melt.h"
#include "material_properties/rho_0.h"
#include "material_properties/P_min.h"

#ifdef MATERIAL_STRENGTH

/**
 * @brief Returns the material parameters of a material.
 *
 * @param mat_id The material ID.
 */
__attribute__((always_inline)) INLINE static const struct mat_params *
material_params(const enum eos_planetary_material_id mat_id) {
  const int mat_index = material_index_from_mat_id(mat_id);
  return &eos.all_mat_params[mat_index];
}

/**
 * @brief Returns the material strength parameters that are loaded from the simulation YAML file.
 */
__attribute__((always_inline)) INLINE static const struct strength_props *
strength_get_props(void) {
  return &eos.strength_props;
}

#ifdef STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG
/**
 * @brief Returns the flaws for tensile damage, loaded from the ICs.
 */
__attribute__((always_inline)) INLINE static const struct strength_flaws *
strength_get_flaws(void) {
  return &eos.strength_flaws;
}
#endif

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

  // Schemes for material properties
  material_T_melt_set_default(mat_params);
  material_rho_0_set_default(mat_params);
  material_P_min_set_default(mat_params);

  // Specific constants for material-strength schemes
  #if defined(STRENGTH_YIELD_STRESS_BENZ_ASPHAUG)
    mat_params->yield_stress_benz_asphaug_Y_0 = 0.f;
  #elif defined(STRENGTH_YIELD_STRESS_COLLINS)
    mat_params->yield_stress_collins_Y_0i = 0.f;
    mat_params->yield_stress_collins_Y_Mi = 0.f;
    mat_params->yield_stress_collins_mu_i = 0.f;
    mat_params->yield_stress_collins_Y_0d = 0.f;
    mat_params->yield_stress_collins_Y_Md = FLT_MAX;
    mat_params->yield_stress_collins_mu_d = 0.f;
  #endif

  #if defined(STRENGTH_YIELD_STRESS_WEAKENING_THERMAL_TANH_MELT)
    mat_params->yield_weakening_thermal_tanh_melt_xi = 1.2f;
  #endif

  #if defined(STRENGTH_YIELD_STRESS_WEAKENING_DENSITY_POWER_LAW)
    mat_params->yield_weakening_density_power_law_mult_param = 0.85f;
    mat_params->yield_weakening_density_power_law_pow_param = 4.f;
  #endif

  #if defined(STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG)
    mat_params->damage_tensile_benz_asphaug_crack_velocity_factor = 0.4f;
  #endif

  #if defined(STRENGTH_DAMAGE_SHEAR_COLLINS)
    mat_params->damage_shear_collins_brittle_to_ductile_pressure = 0.f;
    mat_params->damage_shear_collins_brittle_to_plastic_pressure = 0.f;
    mat_params->damage_shear_collins_plastic_strain_at_failure[0] = 0.01f;
    mat_params->damage_shear_collins_plastic_strain_at_failure[1] = 0.05f;
    mat_params->damage_shear_collins_plastic_strain_at_failure[2] = 0.1f;
    mat_params->damage_shear_collins_plastic_strain_at_failure[3] = 1.f;
  #endif

  #if defined(STRENGTH_ARTIFICIAL_STRESS_BASIS_INDEPENDENT)
    mat_params->artif_stress_basis_independent_transition_factor = 5.f;
  #elif defined(STRENGTH_ARTIFICIAL_STRESS_CARTESIAN_COMPONENTS)
    mat_params->artif_stress_cartesian_components_n = 4.f;
    mat_params->artif_stress_cartesian_components_epsilon = 0.2f;
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

  // Method parameters are now in one section per method
  if ((parser_get_section_id(file_params, "YieldStress") >= 0) ||
      (parser_get_section_id(file_params, "ArtificialStress") >= 0))
    error("%s: the YieldStress and ArtificialStress sections have been split "
          "into one section per method, e.g. YieldStressCollins", param_file);

  // General properties
  mat_params->state_type = (enum mat_state_type)parser_get_opt_param_int(
      file_params, "Material:state_type", mat_params->state_type);
  mat_params->shear_mod = parser_get_opt_param_float(
      file_params, "Material:shear_mod", mat_params->shear_mod);
  mat_params->bulk_mod = parser_get_opt_param_float(
      file_params, "Material:bulk_mod", mat_params->bulk_mod);

  // Schemes for material properties
  material_T_melt_read(mat_params, file_params, param_file);
  material_rho_0_read(mat_params, file_params, param_file);
  material_P_min_read(mat_params, file_params, param_file);

  // Specific constants for material-strength schemes
  #if defined(STRENGTH_YIELD_STRESS_BENZ_ASPHAUG)
    mat_params->yield_stress_benz_asphaug_Y_0 = parser_get_opt_param_float(
        file_params, "YieldStressBenzAsphaug:Y_0", mat_params->yield_stress_benz_asphaug_Y_0);
  #elif defined(STRENGTH_YIELD_STRESS_COLLINS)
    mat_params->yield_stress_collins_Y_0i = parser_get_opt_param_float(
        file_params, "YieldStressCollins:Y_0i", mat_params->yield_stress_collins_Y_0i);
    mat_params->yield_stress_collins_Y_Mi = parser_get_opt_param_float(
        file_params, "YieldStressCollins:Y_Mi", mat_params->yield_stress_collins_Y_Mi);
    mat_params->yield_stress_collins_mu_i = parser_get_opt_param_float(
        file_params, "YieldStressCollins:mu_i", mat_params->yield_stress_collins_mu_i);
    mat_params->yield_stress_collins_Y_0d = parser_get_opt_param_float(
        file_params, "YieldStressCollins:Y_0d", mat_params->yield_stress_collins_Y_0d);
    mat_params->yield_stress_collins_Y_Md = parser_get_opt_param_float(
        file_params, "YieldStressCollins:Y_Md", mat_params->yield_stress_collins_Y_Md);
    mat_params->yield_stress_collins_mu_d = parser_get_opt_param_float(
        file_params, "YieldStressCollins:mu_d", mat_params->yield_stress_collins_mu_d);
  #endif

  #if defined(STRENGTH_YIELD_STRESS_WEAKENING_THERMAL_TANH_MELT)
    mat_params->yield_weakening_thermal_tanh_melt_xi = parser_get_opt_param_float(
        file_params, "YieldStressWeakeningThermalTanhMelt:xi",
        mat_params->yield_weakening_thermal_tanh_melt_xi);
  #endif

  #if defined(STRENGTH_YIELD_STRESS_WEAKENING_DENSITY_POWER_LAW)
    mat_params->yield_weakening_density_power_law_mult_param = parser_get_opt_param_float(
        file_params, "YieldStressWeakeningDensityPowerLaw:mult_param",
        mat_params->yield_weakening_density_power_law_mult_param);
    mat_params->yield_weakening_density_power_law_pow_param = parser_get_opt_param_float(
        file_params, "YieldStressWeakeningDensityPowerLaw:pow_param",
        mat_params->yield_weakening_density_power_law_pow_param);
  #endif

  #if defined(STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG)
    mat_params->damage_tensile_benz_asphaug_crack_velocity_factor = parser_get_opt_param_float(
        file_params, "DamageTensileBenzAsphaug:crack_velocity_factor",
        mat_params->damage_tensile_benz_asphaug_crack_velocity_factor);
  #endif

  #if defined(STRENGTH_DAMAGE_SHEAR_COLLINS)
    mat_params->damage_shear_collins_brittle_to_ductile_pressure = parser_get_opt_param_float(
        file_params, "DamageShearCollins:brittle_to_ductile_pressure",
        mat_params->damage_shear_collins_brittle_to_ductile_pressure);
    mat_params->damage_shear_collins_brittle_to_plastic_pressure = parser_get_opt_param_float(
        file_params, "DamageShearCollins:brittle_to_plastic_pressure",
        mat_params->damage_shear_collins_brittle_to_plastic_pressure);
    parser_get_opt_param_float_array(
        file_params, "DamageShearCollins:plastic_strain_at_failure", 4,
        mat_params->damage_shear_collins_plastic_strain_at_failure);
  #endif

  #if defined(STRENGTH_ARTIFICIAL_STRESS_BASIS_INDEPENDENT)
    mat_params->artif_stress_basis_independent_transition_factor = parser_get_opt_param_float(
        file_params, "ArtificialStressBasisIndependent:transition_factor",
        mat_params->artif_stress_basis_independent_transition_factor);
  #elif defined(STRENGTH_ARTIFICIAL_STRESS_CARTESIAN_COMPONENTS)
    mat_params->artif_stress_cartesian_components_n = parser_get_opt_param_float(
        file_params, "ArtificialStressCartesianComponents:n", mat_params->artif_stress_cartesian_components_n);
    mat_params->artif_stress_cartesian_components_epsilon = parser_get_opt_param_float(
        file_params, "ArtificialStressCartesianComponents:epsilon", mat_params->artif_stress_cartesian_components_epsilon);
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

  // cgs to internal
  mat_params->shear_mod /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
  mat_params->bulk_mod /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);

  // Schemes for material properties
  material_T_melt_convert_units(mat_params, us);
  material_rho_0_convert_units(mat_params, us);
  material_P_min_convert_units(mat_params, us);

  // Specific constants for material-strength schemes
  #if defined(STRENGTH_YIELD_STRESS_BENZ_ASPHAUG)
    // SI to cgs
    mat_params->yield_stress_benz_asphaug_Y_0 *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);

    // cgs to internal
    mat_params->yield_stress_benz_asphaug_Y_0 /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
  #elif defined(STRENGTH_YIELD_STRESS_COLLINS)
    // SI to cgs
    mat_params->yield_stress_collins_Y_0i *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);
    mat_params->yield_stress_collins_Y_Mi *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);
    mat_params->yield_stress_collins_Y_0d *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);
    mat_params->yield_stress_collins_Y_Md *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);

    // cgs to internal
    mat_params->yield_stress_collins_Y_0i /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
    mat_params->yield_stress_collins_Y_Mi /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
    mat_params->yield_stress_collins_Y_0d /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
    mat_params->yield_stress_collins_Y_Md /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
  #endif

  #if defined(STRENGTH_DAMAGE_SHEAR_COLLINS)
    // SI to cgs
    mat_params->damage_shear_collins_brittle_to_ductile_pressure *=
        units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);
    mat_params->damage_shear_collins_brittle_to_plastic_pressure *=
        units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);

    // cgs to internal
    mat_params->damage_shear_collins_brittle_to_ductile_pressure /=
        units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
    mat_params->damage_shear_collins_brittle_to_plastic_pressure /=
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

  if (mat_params->shear_mod <= 0.f)
    error("Material %d: a solid or variable state_type needs "
          "Material:shear_mod > 0", mat_id);

  #if defined(STRENGTH_YIELD_STRESS_COLLINS)
    if ((mat_params->yield_stress_collins_Y_Mi < mat_params->yield_stress_collins_Y_0i) ||
        (mat_params->yield_stress_collins_Y_Md < mat_params->yield_stress_collins_Y_0d))
      error("Material %d: Collins yield stress needs YieldStressCollins:Y_Mi >= Y_0i "
            "and Y_Md >= Y_0d", mat_id);
  #endif

  #if defined(STRENGTH_YIELD_STRESS_WEAKENING_THERMAL_TANH_MELT)
    if (mat_params->yield_weakening_thermal_tanh_melt_xi <= 0.f)
      error("Material %d: thermal yield weakening needs "
            "YieldStressWeakeningThermalTanhMelt:xi > 0", mat_id);
  #endif

  #if defined(STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG)
    if (mat_params->bulk_mod <= 0.f)
      error("Material %d: tensile damage needs Material:bulk_mod > 0", mat_id);
    if (mat_params->damage_tensile_benz_asphaug_crack_velocity_factor <= 0.f)
      error("Material %d: tensile damage needs "
            "DamageTensileBenzAsphaug:crack_velocity_factor > 0", mat_id);
  #endif

  #if defined(STRENGTH_DAMAGE_SHEAR_COLLINS)
    if ((mat_params->damage_shear_collins_brittle_to_ductile_pressure <= 0.f) ||
        (mat_params->damage_shear_collins_brittle_to_plastic_pressure <=
         mat_params->damage_shear_collins_brittle_to_ductile_pressure))
      error("Material %d: shear damage needs 0 < "
            "DamageShearCollins:brittle_to_ductile_pressure < "
            "DamageShearCollins:brittle_to_plastic_pressure", mat_id);
    if ((mat_params->damage_shear_collins_plastic_strain_at_failure[0] <= 0.f) ||
        (mat_params->damage_shear_collins_plastic_strain_at_failure[1] <
         mat_params->damage_shear_collins_plastic_strain_at_failure[0]) ||
        (mat_params->damage_shear_collins_plastic_strain_at_failure[2] <
         mat_params->damage_shear_collins_plastic_strain_at_failure[1]) ||
        (mat_params->damage_shear_collins_plastic_strain_at_failure[3] <
         mat_params->damage_shear_collins_plastic_strain_at_failure[2]))
      error("Material %d: shear damage needs "
            "DamageShearCollins:plastic_strain_at_failure to be positive and "
            "non-decreasing", mat_id);
  #endif

  #if defined(STRENGTH_ARTIFICIAL_STRESS_BASIS_INDEPENDENT)
    if (mat_params->artif_stress_basis_independent_transition_factor <= 0.f)
      error("Material %d: ArtificialStressBasisIndependent:transition_factor must be > 0", mat_id);
  #elif defined(STRENGTH_ARTIFICIAL_STRESS_CARTESIAN_COMPONENTS)
    if ((mat_params->artif_stress_cartesian_components_epsilon > 0.f) &&
        (mat_params->artif_stress_cartesian_components_n <= 0.f))
      error("Material %d: ArtificialStressCartesianComponents:n must be > 0", mat_id);
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

#ifdef MATERIAL_STRENGTH
/**
 * @brief Set the material strength parameters that are loaded from the simulation YAML file.
 *
 * @param sp The material strength parameters.
 * @param params The main parameter file.
 */
INLINE static void set_strength_props(struct strength_props *sp,
                                      struct swift_params *params) {

  sp->strain_rate_timestep_factor = parser_get_opt_param_float(
      params, "Strength:strain_rate_timestep_factor", 0.1f);
  sp->damage_tensile_timestep_factor = parser_get_opt_param_float(
      params, "Strength:damage_tensile_timestep_factor", 0.1f);
  sp->damage_shear_timestep_factor = parser_get_opt_param_float(
      params, "Strength:damage_shear_timestep_factor", 0.1f);
  sp->eta_neighbours = parser_get_param_float(params, "SPH:resolution_eta");

  if (sp->strain_rate_timestep_factor <= 0.f)
    error("Strength:strain_rate_timestep_factor must be > 0");
  if (sp->damage_tensile_timestep_factor <= 0.f)
    error("Strength:damage_tensile_timestep_factor must be > 0");
  if (sp->damage_shear_timestep_factor <= 0.f)
    error("Strength:damage_shear_timestep_factor must be > 0");
}
#endif /* MATERIAL_STRENGTH */

#ifdef STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG
/**
 * @brief Load the flaws for tensile damage from the ICs.
 *
 * @param f The flaws.
 * @param params The main parameter file.
 */
INLINE static void strength_flaws_init(struct strength_flaws *f,
                                       struct swift_params *params) {

#if defined(HAVE_HDF5)
  char ic_file[PARSER_MAX_LINE_SIZE];
  parser_get_param_string(params, "InitialConditions:file_name", ic_file);

  /* Open the file and the group with the flaws */
  const hid_t h_file = H5Fopen(ic_file, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (h_file < 0) error("Error opening the ICs '%s' to read the flaws", ic_file);
  const hid_t h_grp = H5Gopen(h_file, "StrengthFlaws", H5P_DEFAULT);
  if (h_grp < 0) error("The ICs '%s' have no StrengthFlaws group", ic_file);

  /* Get the number of flaws */
  const hid_t h_data = H5Dopen(h_grp, "ActivationThresholds", H5P_DEFAULT);
  if (h_data < 0)
    error("The ICs '%s' have no StrengthFlaws/ActivationThresholds", ic_file);
  f->total_number_of_flaws = io_get_number_element_in_dataset(h_data);
  H5Dclose(h_data);

  /* Read the activation thresholds */
  f->activation_thresholds =
      (float *)malloc(f->total_number_of_flaws * sizeof(float));
  if (f->total_number_of_flaws > 0 && f->activation_thresholds == NULL)
    error("Failed to allocate %lld flaws", f->total_number_of_flaws);
  io_read_array_dataset(h_grp, "ActivationThresholds", FLOAT,
                        f->activation_thresholds, f->total_number_of_flaws);

  H5Gclose(h_grp);
  H5Fclose(h_file);
#else
  error("Reading the flaws for tensile damage needs HDF5");
#endif
}
#endif /* STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG */

#endif /* SWIFT_PLANETARY_MATERIAL_PROPS_H */
