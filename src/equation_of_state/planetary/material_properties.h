/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016   Matthieu Schaller (schaller@strw.leidenuniv.nl).
 *               2018   Jacob Kegerreis (jacob.kegerreis@durham.ac.uk).
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
 * For any/all of the planetary EOS. Each EOS type's functions are set in its
 * own header file: `equation_of_state/planetary/<eos_type>.h`.
 * See `eos_planetary_material_id` for the available choices.
 *
 * Not all functions are implemented for all EOS types, so not all can be used
 * with all hydro formulations yet.
 */

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "common_io.h"
#include "eos_setup.h"
#include "eos_utilities.h"
#include "inline.h"
#include "physical_constants.h"
#include "restart.h"
#include "units.h"

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
 * @brief Returns whether or not the material is in a solid state and not fluid.
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
material_phase_state_from_internal_energy(
    float density, float u, enum eos_planetary_material_id mat_id) {

  const enum eos_planetary_type_id type =
      (enum eos_planetary_type_id)(mat_id / eos_type_factor);
  const int unit_id = mat_id % eos_type_factor;

  const int mat_index = material_index_from_mat_id(mat_id);

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_type_idg:
      return idg_phase_state_from_internal_energy(
          density, u, &eos.mat_params[mat_index], &eos.idg[unit_id]);

    /* Tillotson EoS */
    case eos_type_Til:
      return Til_phase_state_from_internal_energy(
          density, u, &eos.mat_params[mat_index], &eos.Til[unit_id]);

    /* Custom user-provided Tillotson EoS */
    case eos_type_Til_custom:
      return Til_phase_state_from_internal_energy(
          density, u, &eos.mat_params[mat_index], &eos.Til_custom[unit_id]);

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_type_HM80:
      return HM80_phase_state_from_internal_energy(
          density, u, &eos.mat_params[mat_index], &eos.HM80[unit_id]);

    /* SESAME EoS */
    case eos_type_SESAME:
      return SESAME_phase_state_from_internal_energy(
          density, u, &eos.mat_params[mat_index], &eos.SESAME[unit_id]);

    /* ANEOS -- using SESAME-style tables */
    case eos_type_ANEOS:
      return SESAME_phase_state_from_internal_energy(
          density, u, &eos.mat_params[mat_index], &eos.ANEOS[unit_id]);

    /*! Linear EoS -- user-provided parameters */
    case eos_type_linear:
      return linear_phase_state_from_internal_energy(
          density, u, &eos.mat_params[mat_index], &eos.linear[unit_id]);

    /*! Generic user-provided custom tables */
    case eos_type_custom:
      return SESAME_phase_state_from_internal_energy(
          density, u, &eos.mat_params[mat_index], &eos.custom[unit_id]);

    default:
      return -1.f;
  }
}

#ifdef MATERIAL_STRENGTH
/** @brief Returns the shear modulus of a material */
__attribute__((always_inline)) INLINE static float material_shear_mod(
    enum eos_planetary_material_id mat_id) {
  const int mat_index = material_index_from_mat_id(mat_id);
  return eos.mat_params[mat_index].shear_mod;
}

/** @brief Returns the bulk modulus of a material */
__attribute__((always_inline)) INLINE static float material_bulk_mod(
    enum eos_planetary_material_id mat_id) {
  const int mat_index = material_index_from_mat_id(mat_id);
  return eos.mat_params[mat_index].bulk_mod;
}

/** @brief Returns the melting temperature of a material */
__attribute__((always_inline)) INLINE static float material_T_melt(
    enum eos_planetary_material_id mat_id) {
  const int mat_index = material_index_from_mat_id(mat_id);
  return eos.mat_params[mat_index].T_melt;
}

/** @brief Returns the rho_0 of a material */
__attribute__((always_inline)) INLINE static float material_rho_0(
    enum eos_planetary_material_id mat_id) {
  const int mat_index = material_index_from_mat_id(mat_id);
  return eos.mat_params[mat_index].rho_0;
}

// #ifdef STRENGTH_YIELD_###
/** @brief Returns the Y_0 of a material */
__attribute__((always_inline)) INLINE static float material_Y_0(
    enum eos_planetary_material_id mat_id) {
  const int mat_index = material_index_from_mat_id(mat_id);
  return eos.mat_params[mat_index].Y_0;
}

/** @brief Returns the Y_M of a material */
__attribute__((always_inline)) INLINE static float material_Y_M(
    enum eos_planetary_material_id mat_id) {
  const int mat_index = material_index_from_mat_id(mat_id);
  return eos.mat_params[mat_index].Y_M;
}

/** @brief Returns the mu_i of a material */
__attribute__((always_inline)) INLINE static float material_mu_i(
    enum eos_planetary_material_id mat_id) {
  const int mat_index = material_index_from_mat_id(mat_id);
  return eos.mat_params[mat_index].mu_i;
}

/** @brief Returns the mu_d of a material */
__attribute__((always_inline)) INLINE static float material_mu_d(
    enum eos_planetary_material_id mat_id) {
  const int mat_index = material_index_from_mat_id(mat_id);
  return eos.mat_params[mat_index].mu_d;
}

/** @brief Returns the yield stress density softening multiplication parameter
 * of a material */
__attribute__((always_inline)) INLINE static float
material_yield_density_soft_mult_param(enum eos_planetary_material_id mat_id) {
  const int mat_index = material_index_from_mat_id(mat_id);
  return eos.mat_params[mat_index].yield_density_soft_mult_param;
}

/** @brief Returns the yield stress density softening exponent parameter of a
 * material */
__attribute__((always_inline)) INLINE static float
material_yield_density_soft_pow_param(enum eos_planetary_material_id mat_id) {
  const int mat_index = material_index_from_mat_id(mat_id);
  return eos.mat_params[mat_index].yield_density_soft_pow_param;
}

/** @brief Returns the yield stress thermal softening parameter of a material */
__attribute__((always_inline)) INLINE static float
material_yield_thermal_soft_xi(enum eos_planetary_material_id mat_id) {
  const int mat_index = material_index_from_mat_id(mat_id);
  return eos.mat_params[mat_index].yield_thermal_soft_xi;
}

/** @brief Returns the brittle to ductile transition pressure of a material */
__attribute__((always_inline)) INLINE static float
material_brittle_to_ductile_transition_pressure(
    enum eos_planetary_material_id mat_id) {
  const int mat_index = material_index_from_mat_id(mat_id);
  return eos.mat_params[mat_index].brittle_to_ductile_transition_pressure;
}

/** @brief Returns the brittle to plastic transition pressure of a material */
__attribute__((always_inline)) INLINE static float
material_brittle_to_plastic_transition_pressure(
    enum eos_planetary_material_id mat_id) {
  const int mat_index = material_index_from_mat_id(mat_id);
  return eos.mat_params[mat_index].brittle_to_plastic_transition_pressure;
}
// #endif /* STRENGTH_YIELD_### */

#endif /* MATERIAL_STRENGTH */

#endif /* SWIFT_PLANETARY_MATERIAL_PROPS_H */
