/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2024 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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
#ifndef SWIFT_PLANETARY_STRENGTH_PART_DEFAULT_H
#define SWIFT_PLANETARY_STRENGTH_PART_DEFAULT_H

/**
 * @file Planetary/strength/default/hydro_part_strength.h
 * @brief Planetary implementation of SPH with default material strength
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"
#include "strength/damage/damage_data.h"
#include "symmetric_matrix.h"

/**
 * @brief Particle strength fields not needed during the SPH loops over neighbours.
 *
 * This structure contains the particle fields that are not used in the
 * density or force loops. Quantities should be used in the kick, drift and
 * potentially ghost tasks only.
 */
struct strength_xpart_data {

  /*! Evolved density at the last full step. */
  float rho_evol_full;

  // Stress tensor
  struct sym_matrix deviatoric_stress_tensor_full;

  #ifdef STRENGTH_DAMAGE
    // Accumulated damage at the last full step
    struct strength_damage damage_full;
  #endif

};

/**
 * @brief Particle srength fields for the SPH particles
 *
 * The density and force substructures are used to contain variables only used
 * within the density and force loops over neighbours. All more permanent
 * variables should be declared in the main part of the part structure,
 */
struct strength_part_data {

  /*! Particle evolved density (primary density for REMIX equations). */
  float rho_evol;

  /*! Time derivative of the evolved density. */
  float drho_dt;

  // Stress tensor
  struct sym_matrix stress_tensor;

  // Principal stresses (Eigenvalues of stress_tensor)
  float principal_stress_eigen[3];

  // Minimum pressure, from the strength state
  float P_min;

  // Deviatoric stress tensor
  struct sym_matrix deviatoric_stress_tensor;

  // Time derivative of deviatoric stress tensor
  struct sym_matrix dS_dt;

  // Strain rate tensor
  struct sym_matrix strain_rate_tensor;

  // Gradient of velocity, calculated using linear-order reproducing kernel.
  float dv_force_loop[3][3];

#ifdef STRENGTH_DAMAGE
  // Accumulated damage
  struct strength_damage damage;

  // Rates of damage accumulation
  struct strength_damage_rates damage_rates;
#endif

#ifdef STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG
  // Number of flaws for tensile damage accumulation
  int number_of_flaws;

  // Index of the particle's first flaw in the flaws loaded from the ICs
  long long first_flaw_index;
#endif

#ifdef STRENGTH_OBJECT_IDS
  // ID of the object that the particle belongs to
  int object_id;
#endif
};

#endif /* SWIFT_PLANETARY_STRENGTH_PART_DEFAULT_H */
