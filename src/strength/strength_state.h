/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2026 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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
#ifndef SWIFT_STRENGTH_STATE_H
#define SWIFT_STRENGTH_STATE_H

/**
 * @file strength/strength_state.h
 * @brief The state of a particle, as passed to the material strength methods,
 * and the damage of a particle.
 */

#include "symmetric_matrix.h"

/**
 * @brief The damage of a particle.
 */
struct strength_damage {

  /*! Total damage */
  float total;

#ifdef STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG
  /*! Damage accumulated due to tension */
  float tensile;
#endif

#ifdef STRENGTH_DAMAGE_SHEAR_COLLINS
  /*! Damage accumulated due to shear */
  float shear;
#endif
};

/**
 * @brief Particle properties used by the strength methods.
 *
 * Filled in by strength_state_init() with either the drifted or the
 * kicked values of a particle, so the strength methods don't need to read
 * from the particle directly.
 */
struct strength_state {

  /* Set from the particle */

  /*! Phase */
  int phase;

  /*! Mass */
  float mass;

  /*! Density */
  float density;

  /*! Specific internal energy */
  float u;

  /*! Pressure */
  float pressure;

  /*! Temperature */
  float temperature;

  /*! Damage */
  struct strength_damage damage;

  /*! Deviatoric stress tensor */
  struct sym_matrix deviatoric_stress_tensor;

  /*! Strain rate tensor */
  struct sym_matrix strain_rate_tensor;

#ifdef STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG
  /*! Number of flaws for tensile damage accumulation */
  int number_of_flaws;

  /*! Activation thresholds of the flaws (points to the particle's array) */
  const float *activation_thresholds;
#endif

  /* Derived from the above */

  /*! Yield stress, including damage and weakening */
  float yield_stress;

  /*! Stress tensor, including the effects of damage */
  struct sym_matrix stress_tensor;

  /*! Whether the deviatoric stress is on or beyond the yield surface */
  int is_above_yield_criterion;
};

#endif /* SWIFT_STRENGTH_STATE_H */
