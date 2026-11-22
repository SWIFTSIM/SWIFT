/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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
#ifndef SWIFT_PLANETARY_EOS_SETUP_H
#define SWIFT_PLANETARY_EOS_SETUP_H

/**
 * @file equation_of_state/planetary/eos_setup.h
 *
 * Set up elements for the planetary equations of state.
 */

/*! Material identifier flags (material_ID = type_ID * type_factor + unit_ID) */
#define eos_type_factor 100

/**
 * @brief Base type IDs for the planetary EoS.
 */
enum eos_planetary_type_id {

  /*! Ideal gas */
  eos_type_idg = 0,

  /*! Tillotson */
  eos_type_Til = 1,

  /*! Hubbard & MacFarlane (1980) Uranus/Neptune */
  eos_type_HM80 = 2,

  /*! SESAME */
  eos_type_SESAME = 3,

  /*! ANEOS */
  eos_type_ANEOS = 4,

  /*! Linear */
  eos_type_linear = 5,

  /*! Custom */
  eos_type_custom = 9,

  /*! Custom Tillotson */
  eos_type_Til_custom = 10,
};

/**
 * @brief Within-type unit IDs for the planetary EoS.
 */
enum eos_planetary_unit_id {

  /* Ideal gas */

  /*! Default
   * (adiabatic index set at configure time by --with-adiabatic-index, default
   *  value is γ=5/3)
   */
  eos_unit_id_idg_def = 0,

  /* Tillotson */

  /*! Tillotson iron */
  eos_unit_id_Til_iron = 0,

  /*! Tillotson granite */
  eos_unit_id_Til_granite = 1,

  /*! Tillotson water */
  eos_unit_id_Til_water = 2,

  /*! Tillotson basalt */
  eos_unit_id_Til_basalt = 3,

  /*! Tillotson ice */
  eos_unit_id_Til_ice = 4,

  /* Hubbard & MacFarlane (1980) Uranus/Neptune */

  /*! Hydrogen-helium atmosphere */
  eos_unit_id_HM80_HHe = 0,

  /*! H20-CH4-NH3 ice mix */
  eos_unit_id_HM80_ice = 1,

  /*! SiO2-MgO-FeS-FeO rock mix */
  eos_unit_id_HM80_rock = 2,

  /* SESAME (and SESAME-style) */

  /*! SESAME iron 2140 */
  eos_unit_id_SESAME_iron = 0,

  /*! SESAME basalt 7530 */
  eos_unit_id_SESAME_basalt = 1,

  /*! SESAME water 7154 */
  eos_unit_id_SESAME_water = 2,

  /*! Senft & Stewart (2008) water */
  eos_unit_id_SS08_water = 3,

  /*! AQUA water (Haldemann et al. 2020) */
  eos_unit_id_AQUA = 4,

  /*! CMS19 (Chabrier et al. 2019) hydrogen */
  eos_unit_id_CMS19_H = 5,

  /*! CMS19 (Chabrier et al. 2019) helium */
  eos_unit_id_CMS19_He = 6,

  /*! CMS19 (Chabrier et al. 2019) H--He mixture (Y=0.245) */
  eos_unit_id_CD21_HHe = 7,

  /* ANEOS */

  /*! ANEOS forsterite (Stewart et al. 2019) -- in SESAME-style tables */
  eos_unit_id_ANEOS_forsterite = 0,

  /*! ANEOS iron (Stewart 2020) -- in SESAME-style tables */
  eos_unit_id_ANEOS_iron = 1,

  /*! ANEOS Fe85Si15 (Stewart 2020) -- in SESAME-style tables */
  eos_unit_id_ANEOS_Fe85Si15 = 2,
};

/**
 * @brief Material IDs (combined type and unit) for the planetary EoS.
 */
enum eos_planetary_material_id {

  /* Ideal gas */

  /*! Default
   * (adiabatic index set at configure time by --with-adiabatic-index, default
   *  value is γ=5/3)
   */
  eos_mat_id_idg_def = eos_type_idg * eos_type_factor + eos_unit_id_idg_def,

  /* Tillotson */

  /*! Tillotson iron */
  eos_mat_id_Til_iron = eos_type_Til * eos_type_factor + eos_unit_id_Til_iron,

  /*! Tillotson granite */
  eos_mat_id_Til_granite =
      eos_type_Til * eos_type_factor + eos_unit_id_Til_granite,

  /*! Tillotson water */
  eos_mat_id_Til_water = eos_type_Til * eos_type_factor + eos_unit_id_Til_water,

  /*! Tillotson basalt */
  eos_mat_id_Til_basalt =
      eos_type_Til * eos_type_factor + eos_unit_id_Til_basalt,

  /*! Tillotson ice */
  eos_mat_id_Til_ice = eos_type_Til * eos_type_factor + eos_unit_id_Til_ice,

  /* Hubbard & MacFarlane (1980) Uranus/Neptune */

  /*! Hydrogen-helium atmosphere */
  eos_mat_id_HM80_HHe = eos_type_HM80 * eos_type_factor + eos_unit_id_HM80_HHe,

  /*! H20-CH4-NH3 ice mix */
  eos_mat_id_HM80_ice = eos_type_HM80 * eos_type_factor + eos_unit_id_HM80_ice,

  /*! SiO2-MgO-FeS-FeO rock mix */
  eos_mat_id_HM80_rock =
      eos_type_HM80 * eos_type_factor + eos_unit_id_HM80_rock,

  /* SESAME (and SESAME-style) */

  /*! SESAME iron 2140 */
  eos_mat_id_SESAME_iron =
      eos_type_SESAME * eos_type_factor + eos_unit_id_SESAME_iron,

  /*! SESAME basalt 7530 */
  eos_mat_id_SESAME_basalt =
      eos_type_SESAME * eos_type_factor + eos_unit_id_SESAME_basalt,

  /*! SESAME water 7154 */
  eos_mat_id_SESAME_water =
      eos_type_SESAME * eos_type_factor + eos_unit_id_SESAME_water,

  /*! Senft & Stewart (2008) water */
  eos_mat_id_SS08_water =
      eos_type_SESAME * eos_type_factor + eos_unit_id_SS08_water,

  /*! AQUA water (Haldemann et al. 2020) */
  eos_mat_id_AQUA = eos_type_SESAME * eos_type_factor + eos_unit_id_AQUA,

  /*! CMS19 (Chabrier et al. 2019) hydrogen */
  eos_mat_id_CMS19_H = eos_type_SESAME * eos_type_factor + eos_unit_id_CMS19_H,

  /*! CMS19 (Chabrier et al. 2019) helium */
  eos_mat_id_CMS19_He =
      eos_type_SESAME * eos_type_factor + eos_unit_id_CMS19_He,

  /*! CMS19 (Chabrier et al. 2019) H--He mixture (Y=0.245) */
  eos_mat_id_CD21_HHe =
      eos_type_SESAME * eos_type_factor + eos_unit_id_CD21_HHe,

  /* ANEOS */

  /*! ANEOS forsterite (Stewart et al. 2019) -- in SESAME-style tables */
  eos_mat_id_ANEOS_forsterite =
      eos_type_ANEOS * eos_type_factor + eos_unit_id_ANEOS_forsterite,

  /*! ANEOS iron (Stewart 2020) -- in SESAME-style tables */
  eos_mat_id_ANEOS_iron =
      eos_type_ANEOS * eos_type_factor + eos_unit_id_ANEOS_iron,

  /*! ANEOS Fe85Si15 (Stewart 2020) -- in SESAME-style tables */
  eos_mat_id_ANEOS_Fe85Si15 =
      eos_type_ANEOS * eos_type_factor + eos_unit_id_ANEOS_Fe85Si15,
};

/**
 * @brief State type of the material.
 */
enum mat_state_type {
  /*! Always fluid */
  mat_state_type_fluid = 0,

  /*! Always solid */
  mat_state_type_solid = 1,

  /*! Variable */
  mat_state_type_variable = 2,
};

/**
 * @brief Phase state of a particle.
 */
enum mat_phase {
  mat_phase_unset = 0,
  mat_phase_fluid = 1,
  mat_phase_solid = 2,
};

#ifdef MATERIAL_STRENGTH
/**
 * @brief Schemes for the melting temperature of a material.
 */
enum T_melt_scheme {
  T_melt_scheme_none = 0,
  T_melt_scheme_constant = 1,
};

/**
 * @brief Schemes for the reference density of a material.
 */
enum rho_0_scheme {
  rho_0_scheme_none = 0,
  rho_0_scheme_constant = 1,
  rho_0_scheme_from_eos = 2,
};

/**
 * @brief Schemes for the minimum pressure of a material.
 */
enum P_min_scheme {
  P_min_scheme_none = 0,
  P_min_scheme_constant = 1,
  P_min_scheme_yield_envelope = 2,
};
#endif /* MATERIAL_STRENGTH */

/**
 * @brief Struct of material parameters beyond the base EoS.
 */
struct mat_params {
  enum mat_state_type state_type;

#ifdef MATERIAL_STRENGTH
  float shear_mod;
  float bulk_mod;
  enum T_melt_scheme T_melt_scheme;
  float T_melt_constant;
  enum rho_0_scheme rho_0_scheme;
  float rho_0_constant;
  enum P_min_scheme P_min_scheme;
  float P_min_constant;

#if defined(STRENGTH_YIELD_STRESS_BENZ_ASPHAUG)
  float yield_stress_benz_asphaug_Y_0;
#elif defined(STRENGTH_YIELD_STRESS_COLLINS)
  float yield_stress_collins_Y_0i;
  float yield_stress_collins_Y_Mi;
  float yield_stress_collins_mu_i;
  float yield_stress_collins_Y_0d;
  float yield_stress_collins_Y_Md;
  float yield_stress_collins_mu_d;
#endif

#if defined(STRENGTH_YIELD_STRESS_WEAKENING_THERMAL_TANH_MELT)
  float yield_weakening_thermal_tanh_melt_xi;
#endif

#if defined(STRENGTH_YIELD_STRESS_WEAKENING_DENSITY_POWER_LAW)
  float yield_weakening_density_power_law_mult_param;
  float yield_weakening_density_power_law_pow_param;
#endif

#if defined(STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG)
  float damage_tensile_benz_asphaug_crack_velocity_factor;
#endif

#if defined(STRENGTH_DAMAGE_SHEAR_COLLINS)
  float damage_shear_collins_brittle_to_ductile_pressure;
  float damage_shear_collins_brittle_to_plastic_pressure;
  float damage_shear_collins_plastic_strain_at_failure[4];
#endif

#if defined(STRENGTH_ARTIFICIAL_STRESS_BASIS_INDEPENDENT)
  float artif_stress_basis_independent_transition_factor;
#elif defined(STRENGTH_ARTIFICIAL_STRESS_CARTESIAN_COMPONENTS)
  float artif_stress_cartesian_components_n;
  float artif_stress_cartesian_components_epsilon;
#endif
#endif /* MATERIAL_STRENGTH */
};

#ifdef MATERIAL_STRENGTH
/**
 * @brief Material strength parameters that are loaded from the simulation YAML
 * file.
 */
struct strength_props {

  /*! Maximum growth in cbrt(tensile damage) per time-step */
  float damage_tensile_timestep_factor;

  /*! Maximum growth in shear damage per time-step */
  float damage_shear_timestep_factor;

  /*! Target smoothing length in units of the mean inter-particle separation */
  float eta_neighbours;
};
#endif /* MATERIAL_STRENGTH */

#ifdef STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG
/**
 * @brief Flaws for tensile damage, loaded from the ICs.
 *
 * The activation thresholds of every particle's flaws, one particle after
 * another. Each particle stores the index of its first flaw and its number of
 * flaws.
 */
struct strength_flaws {

  /*! Total number of flaws */
  long long total_number_of_flaws;

  /*! Activation thresholds of all the flaws */
  float *activation_thresholds;
};
#endif /* STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG */

/**
 * @brief Count the number of materials of each EoS type.
 */
enum eos_type_count {
  eos_count_idg = 1,
  eos_count_Til = 5,
  eos_count_Til_custom = 10,
  eos_count_HM80 = 3,
  eos_count_SESAME = 8,
  eos_count_ANEOS = 3,
  eos_count_linear = 10,
  eos_count_custom = 10,
};

/**
 * @brief Cumulative consecutive count (starting index) for each material type.
 */
enum eos_type_cumul_count {
  eos_cumul_count_idg = 0,
  eos_cumul_count_Til = eos_cumul_count_idg + eos_count_idg,
  eos_cumul_count_Til_custom = eos_cumul_count_Til + eos_count_Til,
  eos_cumul_count_HM80 = eos_cumul_count_Til_custom + eos_count_Til_custom,
  eos_cumul_count_SESAME = eos_cumul_count_HM80 + eos_count_HM80,
  eos_cumul_count_ANEOS = eos_cumul_count_SESAME + eos_count_SESAME,
  eos_cumul_count_linear = eos_cumul_count_ANEOS + eos_count_ANEOS,
  eos_cumul_count_custom = eos_cumul_count_linear + eos_count_linear,

  // Total count of materials
  eos_count_total = eos_cumul_count_custom + eos_count_custom + 1,
};

#endif /* SWIFT_PLANETARY_EOS_SETUP_H */
