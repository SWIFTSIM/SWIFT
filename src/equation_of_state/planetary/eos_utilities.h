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
#include "eos_setup.h"
#include "eos_hm80.h"
#include "eos_ideal_gas.h"
#include "eos_linear.h"
#include "eos_sesame.h"
#include "eos_tillotson.h"

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

/*
    Read the extra material parameters from a file. ###change to .yml or
   something

    File contents
    -------------
    # header (2 lines)
    phase_state (enum mat_phase_state: fluid=0, solid=1, variable=2)  shear_mod
   (Pa)
*/
INLINE static void set_material_params(struct mat_params *all_mat_params,
                                       enum eos_planetary_material_id mat_id,
                                       char *param_file,
                                       const struct unit_system *us) {

  const int mat_index = material_index_from_mat_id(mat_id);

#ifdef MATERIAL_STRENGTH
  // Load table contents from file
  FILE *f = fopen(param_file, "r");
  if (f == NULL)
    error("Failed to open the material parameters file '%s'", param_file);

  // Skip header lines
  skip_lines(f, 2);

  // Read parameters
  int c, phase_state;
  c = fscanf(f, "%d %f", &phase_state, &all_mat_params[mat_index].shear_mod);
  if (c != 2)
    error("Failed to read the material parameters file %s", param_file);

  all_mat_params[mat_index].phase_state = (enum mat_phase_state)phase_state;
#else
  all_mat_params[mat_index].phase_state = mat_phase_state_fluid;
#endif /* MATERIAL_STRENGTH */

  // ###convert units!
  //  us...
}

#endif /* SWIFT_PLANETARY_EOS_UTILITIES_H */
