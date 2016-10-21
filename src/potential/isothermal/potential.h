/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Tom Theuns (tom.theuns@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_POTENTIAL_ISOTHERMAL_H
#define SWIFT_POTENTIAL_ISOTHERMAL_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <math.h>

/* Local includes. */
#include "error.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"

/**
 * @brief External Potential Properties - Isothermal sphere case
 */
struct external_potential {

  /*! Position of the centre of potential */
  double x, y, z;

  /*! Rotation velocity */
  double vrot;

  /*! Square of vrot divided by G \f$ \frac{v_{rot}^2}{G} \f$ */
  double vrot2_over_G;

  /*! Time-step condition pre-factor */
  double timestep_mult;
};

/**
 * @brief Computes the time-step due to the acceleration from an isothermal
 * potential.
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

  const float dx = g->x[0] - potential->x;
  const float dy = g->x[1] - potential->y;
  const float dz = g->x[2] - potential->z;

  const float rinv2 = 1.f / (dx * dx + dy * dy + dz * dz);
  const float drdv =
      dx * (g->v_full[0]) + dy * (g->v_full[1]) + dz * (g->v_full[2]);
  const double vrot = potential->vrot;

  const float dota_x =
      vrot * vrot * rinv2 * (g->v_full[0] - 2.f * drdv * dx * rinv2);
  const float dota_y =
      vrot * vrot * rinv2 * (g->v_full[1] - 2.f * drdv * dy * rinv2);
  const float dota_z =
      vrot * vrot * rinv2 * (g->v_full[2] - 2.f * drdv * dz * rinv2);
  const float dota_2 = dota_x * dota_x + dota_y * dota_y + dota_z * dota_z;
  const float a_2 = g->a_grav[0] * g->a_grav[0] + g->a_grav[1] * g->a_grav[1] +
                    g->a_grav[2] * g->a_grav[2];

  return potential->timestep_mult * sqrtf(a_2 / dota_2);
}

/**
 * @brief Computes the gravitational acceleration from an isothermal potential.
 *
 * Note that the accelerations are multiplied by Newton's G constant
 * later on.
 *
 * @param time The current time.
 * @param potential The #external_potential used in the run.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static void external_gravity_acceleration(
    double time, const struct external_potential* potential,
    const struct phys_const* const phys_const, struct gpart* g) {

  const float dx = g->x[0] - potential->x;
  const float dy = g->x[1] - potential->y;
  const float dz = g->x[2] - potential->z;

  const float rinv2 = 1.f / (dx * dx + dy * dy + dz * dz);

  const double term = -potential->vrot2_over_G * rinv2;

  g->a_grav[0] += term * dx;
  g->a_grav[1] += term * dy;
  g->a_grav[2] += term * dz;
}

/**
 * @brief Computes the gravitational potential energy of a particle in an
 * isothermal potential.
 *
 * @param potential The #external_potential used in the run.
 * @param phys_const Physical constants in internal units.
 * @param g Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
external_gravity_get_potential_energy(
    const struct external_potential* potential,
    const struct phys_const* const phys_const, const struct gpart* g) {

  const float dx = g->x[0] - potential->x;
  const float dy = g->x[1] - potential->y;
  const float dz = g->x[2] - potential->z;

  return 0.5f * potential->vrot * potential->vrot *
         logf(dx * dx + dy * dy * dz * dz);
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
    const struct swift_params* parameter_file,
    const struct phys_const* phys_const, const struct UnitSystem* us,
    const struct space* s, struct external_potential* potential) {

  potential->x =
      s->dim[0] / 2. +
      parser_get_param_double(parameter_file, "IsothermalPotential:position_x");
  potential->y =
      s->dim[1] / 2. +
      parser_get_param_double(parameter_file, "IsothermalPotential:position_y");
  potential->z =
      s->dim[2] / 2. +
      parser_get_param_double(parameter_file, "IsothermalPotential:position_z");
  potential->vrot =
      parser_get_param_double(parameter_file, "IsothermalPotential:vrot");
  potential->timestep_mult = parser_get_param_float(
      parameter_file, "IsothermalPotential:timestep_mult");

  potential->vrot2_over_G =
      potential->vrot * potential->vrot / phys_const->const_newton_G;
}

/**
 * @brief Prints the properties of the external potential to stdout.
 *
 * @param  potential The external potential properties.
 */
static INLINE void potential_print_backend(
    const struct external_potential* potential) {

  message(
      "External potential is 'Isothermal' with properties are (x,y,z) = (%e, "
      "%e, %e), vrot = %e "
      "timestep multiplier = %e.",
      potential->x, potential->y, potential->z, potential->vrot,
      potential->timestep_mult);
}

#endif /* SWIFT_POTENTIAL_ISOTHERMAL_H */
