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
#ifndef SWIFT_REMIX_STRENGTH_PART_DEFAULT_H
#define SWIFT_REMIX_STRENGTH_PART_DEFAULT_H

/**
 * @file REMIX/strength/default/hydro_part_strength.h
 * @brief REMIX implementation of SPH with default material strength
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"
#include "symmetric_matrix.h"

/**
 * @brief Particle strength fields not needed during the SPH loops over neighbours.
 *
 * This structure contains the particle fields that are not used in the
 * density or force loops. Quantities should be used in the kick, drift and
 * potentially ghost tasks only.
 */
struct strength_xpart_data {  

  // Stress tensor
  struct sym_matrix deviatoric_stress_tensor_full;

#if defined(STRENGTH_DAMAGE)
  // Accumulated damage at the last full step
  float damage_full;

  // Damage accumulated due to tension at the last full step
  float tensile_damage_full;

  // Damage accumulated due to shear at the last full step
  float shear_damage_full;
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

  // Stress tensor
  struct sym_matrix stress_tensor;

  // Principal stresses (Eigenvalues of stress_tensor)
  float principal_stress_eigen[3];

  // Deviatoric stress tensor
  struct sym_matrix deviatoric_stress_tensor;

  // Time derivative of deviatoric stress tensor
  struct sym_matrix dS_dt;

  // Gradient of velocity, calculated using linear-order reproducing kernel.
  float dv_force_loop[3][3];

#if defined(STRENGTH_DAMAGE)
  // Accumulated damage
  float damage;

  // Damage accumulated due to tension
  float tensile_damage;

  // Damage accumulated due to shear
  float shear_damage;

  // Need to store this as a particle parameter for timestep
  float dD_dt;
#endif

#if defined(STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG)
  // Number of flaws for tensile damage accumulation
  int number_of_flaws;

  // Activation thresholds of flaws for tensile damage accumulation
  // ### Work out how to set the length of this
  float activation_thresholds[40];
#endif
};

#endif /* SWIFT_REMIX_STRENGTH_PART_DEFAULT_H */