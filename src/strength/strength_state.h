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

#include <stdbool.h>

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
 * @brief Quantities derived from the strength state.
 *
 * Each one is computed the first time it is asked for, by its
 * strength_state_get_*() function, and then reused.
 */
struct strength_state_cache {

  /*! Pressure */
  bool pressure_computed;
  float pressure;

  /*! Temperature */
  bool temperature_computed;
  float temperature;

  /*! J_2 invariant of the deviatoric stress tensor */
  bool J_2_computed;
  float J_2;

  /*! Stress tensor, including the effects of damage */
  bool stress_tensor_computed;
  struct sym_matrix stress_tensor;

#ifdef STRENGTH_YIELD_STRESS
  /*! Yield stress, including damage and weakening */
  bool yield_stress_computed;
  float yield_stress;
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

  /*! Material ID */
  enum eos_planetary_material_id mat_id;

  /*! Phase */
  int phase;

  /*! Mass */
  float mass;

  /*! Density */
  float density;

  /*! Specific internal energy */
  float u;

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

  /*! Quantities derived from the above, only read through
   * strength_state_get_*() */
  struct strength_state_cache cache;
};

/*
 * Getters for the derived quantities. They are declared here so that any
 * strength method can use them, and defined in strength_state_setup.h or, for
 * quantities that depend on a method, in that method's file.
 */
__attribute__((always_inline)) INLINE static float strength_state_get_pressure(
    struct strength_state *s, const struct mat_params *mp);

__attribute__((always_inline)) INLINE static float strength_state_get_temperature(
    struct strength_state *s, const struct mat_params *mp);

__attribute__((always_inline)) INLINE static float strength_state_get_J_2(
    struct strength_state *s, const struct mat_params *mp);

__attribute__((always_inline)) INLINE static struct sym_matrix
strength_state_get_stress_tensor(struct strength_state *s,
                                 const struct mat_params *mp);

#ifdef STRENGTH_YIELD_STRESS
__attribute__((always_inline)) INLINE static float strength_state_get_yield_stress(
    struct strength_state *s, const struct mat_params *mp);
#endif

#endif /* SWIFT_STRENGTH_STATE_H */
