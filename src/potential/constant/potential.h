/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019  Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_POTENTIAL_CONSTANT_H
#define SWIFT_POTENTIAL_CONSTANT_H

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <float.h>
#include <math.h>

/* Local includes. */
#include "error.h"
#include "gravity.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"

/**
 * @brief External Potential Properties - Constant acceleration.
 */
struct external_potential {

  /*! Value of the acceleration */
  double g[3];
};

/**
 * @brief Computes the time-step due to the acceleration from an constant
 * acceleration.
 *
 * @param time The current time.
 * @param potential The #external_potential used in the run.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static float external_gravity_timestep(
    double time, const struct external_potential* restrict potential,
    const struct phys_const* restrict phys_const,
    const struct gpart* restrict g) {

  return FLT_MAX;
}

/**
 * @brief Computes the gravitational acceleration from an constant acceleration.
 *
 * @param time The current time.
 * @param potential The #external_potential used in the run.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static void external_gravity_acceleration(
    double time, const struct external_potential* potential,
    const struct phys_const* const phys_const, struct gpart* g) {

  const float gh = g->x[0] * potential->g[0] + g->x[1] * potential->g[1] +
                   g->x[2] * potential->g[2];
  const float pot = -gh / 3.f;

  g->a_grav[0] += potential->g[0];
  g->a_grav[1] += potential->g[1];
  g->a_grav[2] += potential->g[2];
  gravity_add_comoving_potential(g, pot);
}

/**
 * @brief Computes the gravitational potential energy of a particle in an
 * constant acceleration.
 *
 * @param time The current time (unused here).
 * @param potential The #external_potential used in the run.
 * @param phys_const Physical constants in internal units.
 * @param g Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
external_gravity_get_potential_energy(
    double time, const struct external_potential* potential,
    const struct phys_const* const phys_const, const struct gpart* g) {

  const float gh = g->x[0] * potential->g[0] + g->x[1] * potential->g[1] +
                   g->x[2] * potential->g[2];

  return g->mass * gh;
}

/**
 * @brief Initialises the external potential properties in the internal system
 * of units.
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param potential The external potential properties to initialize
 */
static INLINE void potential_init_backend(
    struct swift_params* parameter_file, const struct phys_const* phys_const,
    const struct unit_system* us, const struct space* s,
    struct external_potential* potential) {

  /* Read in the acceleration */
  parser_get_param_double_array(parameter_file, "ConstantPotential:g_cgs", 3,
                                potential->g);

  /* Change the unit system */
  const double unit_length = units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
  const double unit_time = units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  const double unit_g = unit_length / (unit_time * unit_time);

  for (int i = 0; i < 3; i++) {
    // Need to divide by G due to gravity_end_force
    potential->g[i] /= unit_g * phys_const->const_newton_G;
  }
}

/**
 * @brief Prints the properties of the external potential to stdout.
 *
 * @param  potential The external potential properties.
 */
static INLINE void potential_print_backend(
    const struct external_potential* potential) {

  message(
      "External potential is 'Constant' with properties are g = (%e, "
      "%e, %e)",
      potential->g[0], potential->g[1], potential->g[2]);
}

#endif /* SWIFT_CONSTANT_H */
