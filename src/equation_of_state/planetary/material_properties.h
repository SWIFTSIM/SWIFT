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
#include "inline.h"
#include "physical_constants.h"
#include "restart.h"
#include "units.h"

#include "equation_of_state.h"

/**
 * @brief Phase state of the material.
 */
enum mat_phase_state {
  /*! Always fluid */
  mat_phase_state_fluid = 0,

  /*! Always solid */
  mat_phase_state_solid = 1,

  /*! Variable */
  mat_phase_state_variable = 2,
};

/**
 * @brief Material parameters beyond the base EoS
 */
struct mat_params {
  enum eos_phase_state phase_state;

#ifdef MATERIAL_STRENGTH
  float shear_mod;
  float bulk_mod;
  float T_melt;

#if defined(STRENGTH_YIELD_###)
  float Y_0;
  float Y_M;
#endif

#endif /* MATERIAL_STRENGTH */
};

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
      return mat_index_base_idg + unit_id;

    /* Tillotson EoS */
    case eos_type_Til:
      return mat_index_base_Til + unit_id;

    /* Custom user-provided Tillotson EoS */
    case eos_type_Til_custom:
      return mat_index_base_Til_custom + unit_id;

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_type_HM80:
      return mat_index_base_HM80 + unit_id;

    /* SESAME EoS */
    case eos_type_SESAME:
      return mat_index_base_SESAME + unit_id;

    /* ANEOS -- using SESAME-style tables */
    case eos_type_ANEOS:
      return mat_index_base_ANEOS + unit_id;

    /*! Linear EoS -- user-provided parameters */
    case eos_type_linear:
      return mat_index_base_linear + unit_id;

    /*! Generic user-provided custom tables */
    case eos_type_custom:
      return mat_index_base_custom + unit_id;

    default:
      return -1.f;
  }
}

/*
    Read the extra material parameters from a file. ###change to .yml or something

    File contents
    -------------
    # header (2 lines)
    phase_state (enum eos_phase_state: fluid=0, solid=1, variable=2)  shear_mod (Pa)
*/
INLINE static void set_material_params(
    struct material_parameters *mat_params, enum eos_planetary_material_id mat_id,
    char *param_file, const struct unit_system *us) {

  const int mat_index = material_index_from_mat_id(mat_id);

#ifdef MATERIAL_STRENGTH
  // Load table contents from file
  FILE *f = fopen(param_file, "r");
  if (f == NULL)
    error("Failed to open the material parameters file '%s'", param_file);

  // Skip header lines
  skip_lines(f, 2);

  // Read parameters
  int c;
  c = fscanf(f, "%d %f", &mat_params[i_mat].phase_state, &mat_params[i_mat].shear_mod);
  if (c != 2) error("Failed to read the material parameters file %s", param_file);
#else
  mat_params[i_mat].phase_state = eos_phase_state_fluid;
#endif /* MATERIAL_STRENGTH */

  //###convert units!
  // us...
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

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_type_idg:
      return idg_phase_state_from_internal_energy(density, u, &eos.idg[unit_id]);

    /* Tillotson EoS */
    case eos_type_Til:
      return Til_phase_state_from_internal_energy(density, u, &eos.Til[unit_id]);

    /* Custom user-provided Tillotson EoS */
    case eos_type_Til_custom:
      return Til_phase_state_from_internal_energy(density, u, &eos.Til_custom[unit_id]);

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_type_HM80:
      return HM80_phase_state_from_internal_energy(density, u, &eos.HM80[unit_id]);

    /* SESAME EoS */
    case eos_type_SESAME:
      return SESAME_phase_state_from_internal_energy(density, u, &eos.SESAME[unit_id]);

    /* ANEOS -- using SESAME-style tables */
    case eos_type_ANEOS:
      return SESAME_phase_state_from_internal_energy(density, u, &eos.ANEOS[unit_id]);

    /*! Linear EoS -- user-provided parameters */
    case eos_type_linear:
      return linear_phase_state_from_internal_energy(density, u, &eos.linear[unit_id]);

    /*! Generic user-provided custom tables */
    case eos_type_custom:
      return SESAME_phase_state_from_internal_energy(density, u, &eos.custom[unit_id]);

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
#endif /* MATERIAL_STRENGTH */

#ifdef STRENGTH_YIELD_###
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
#endif /* STRENGTH_YIELD_### */

#endif /* SWIFT_PLANETARY_MATERIAL_PROPS_H */
