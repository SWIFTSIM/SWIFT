/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_ADAPTIVE_SOFTENING_STRUCT_H
#define SWIFT_ADAPTIVE_SOFTENING_STRUCT_H

/* Config parameters. */
#include <config.h>

#ifdef ADAPTIVE_SOFTENING

/**
 * @brief Particle-carried fields for the adaptive softening scheme.
 */
struct adaptive_softening_part_data {

  /*! Correction term for energy conservation */
  float zeta;
};

/**
 * @brief Gravity particle-carried fields for the adaptive softening scheme.
 */

/* 1) Compute zeta and omega during the tree walk with the tidal tensor of the
  previous timestep. We can drift them to et an estimate of the current
  timestep, but we'll do this later. At the same time, compute the new tidal
  tensor. (store the old value in another variable.
  ---> create a new adaptive_softening_add_correction_term() fct and call it
       in runner_iact_grav_pp_full() and _truncated().

       To finish the computations, we might need to add a task.

  2) Compute the Gamma_ab terms during the force interaction, i.e. when we
  compute the acceleration. Compute d_epsilon_dt. Pay attention to the
  interaction with gas particles -> not the same softening rule.
  --> create a new adaptive_softening_get_acc_term() and call it in

  3) Update the softening with the new tidal tensor. We can do this probably in
  grav_end_force or in kick2() ? See where the gas updates its value.


  Notes:
  We only need to compute zeta and omega over all particles in the softening
  kernel. Outside the softening kernel (r>H), the d Phi/d_epsilon = 0.

  When computing the tidal tensor, use do not use the same Phi. See Hopkins
  paper.

  At the end, add a minimal softening and maximal softening.
*/
struct adaptive_softening_gpart_data {

  /*! Correction term for energy conservation eq 15 in Hopkins 2023.*/
  float zeta;

  /*! Correction term for energy conservation (second term) eq 35 */
  float omega;

  /*! Time derivative of the softening eq 37 */
  float depsilon_dt;
};

#else

/**
 * @brief Particle-carried fields for the adaptive softening scheme.
 */
struct adaptive_softening_part_data {};

/**
 * @brief Gravity particle-carried fields for the adaptive softening scheme.
 */
struct adaptive_softening_gpart_data {};

#endif

#endif /* SWIFT_ADAPTIVE_SOFTENING_STRUCT_H */
