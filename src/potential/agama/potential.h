/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Tom Theuns (tom.theuns@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_POTENTIAL_AGAMA_H
#define SWIFT_POTENTIAL_AGAMA_H

/* Config parameters. */
#include <config.h>
#include <agama/interface_c.h>

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
 * @brief External Potential Properties - Point mass case
 */
struct external_potential {

  /* Potential object */
  agama_Potential *agama_potential;

  /*! Position of the point mass */
  double x[3];

  /*! Time-step condition pre-factor */
  float timestep_mult;
};

/**
 * @brief Computes the time-step due to the acceleration from a point mass
 *
 * We pass in the time for simulations where the potential evolves with time.
 *
 * @param time The current time.
 * @param potential The properties of the externa potential.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static float external_gravity_timestep(
    double time, const struct external_potential *restrict potential,
    const struct phys_const *restrict phys_const,
    const struct gpart *restrict g) {

  const float a_2 = g->a_grav[0]*g->a_grav[0] + g->a_grav[1]*g->a_grav[1] + g->a_grav[2]*g->a_grav[2];
  const float v_2 = g->v_full[0]*g->v_full[0] + g->v_full[1]*g->v_full[1] + g->v_full[2]*g->v_full[2];

  if (fabsf(a_2) > 0.f)
    return potential->timestep_mult * sqrtf(v_2 / a_2);
  else
    return FLT_MAX;
}

/**
 * @brief Computes the gravitational acceleration of a particle due to a
 * point mass
 *
 * Note that the accelerations are multiplied by Newton's G constant later
 * on.
 *
 * We pass in the time for simulations where the potential evolves with time.
 *
 * @param time The current time.
 * @param potential The proerties of the external potential.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static void external_gravity_acceleration(
    double time, const struct external_potential *restrict potential,
    const struct phys_const *restrict phys_const, struct gpart *restrict g) {


  double pot;
  double t=0;
  double deriv[3];
  double deriv2[6];

  /* Center the potential */
  const float dx = g->x[0] - potential->x[0];
  const float dy = g->x[1] - potential->x[1];
  const float dz = g->x[2] - potential->x[2];
  double xyz[3] = {dx, dy, dz};

  /* Compute potential and forces */  
  pot = agama_evalPotential(potential->agama_potential, xyz, t, deriv, deriv2);

  g->a_grav[0] += -deriv[0];
  g->a_grav[1] += -deriv[1];
  g->a_grav[2] += -deriv[2];
  gravity_add_comoving_potential(g, pot);
}

/**
 * @brief Computes the gravitational potential energy of a particle in a point
 * mass potential.
 *
 * @param time The current time (unused here).
 * @param potential The #external_potential used in the run.
 * @param phys_const Physical constants in internal units.
 * @param g Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
external_gravity_get_potential_energy(
    double time, const struct external_potential *potential,
    const struct phys_const *const phys_const, const struct gpart *g) {

  double pot;
  double t=0;
  double deriv[3];
  double deriv2[6];
  
  /* Center the potential */
  const float dx = g->x[0] - potential->x[0];
  const float dy = g->x[1] - potential->x[1];
  const float dz = g->x[2] - potential->x[2];
  double xyz[3] = {dx, dy, dz};

  /* Compute potential and forces */  
  pot = agama_evalPotential(potential->agama_potential, xyz, t, deriv, deriv2);
  
  return pot;
}

/**
 * @brief Initialises the external potential properties in the internal system
 * of units.
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param s The #space we run in.
 * @param potential The external potential properties to initialize
 */
static INLINE void potential_init_backend(
    struct swift_params *parameter_file, const struct phys_const *phys_const,
    const struct unit_system *us, const struct space *s,
    struct external_potential *potential) {



  /* Read in the name of the AGAMA potential file */
  char agama_potential_filename[32];
  parser_get_param_string(parameter_file, "AGAMAPotential:filename",
                                agama_potential_filename);


  /* Read in the position of the centre of potential */
  parser_get_param_double_array(parameter_file, "AGAMAPotential:position",
                                3, potential->x);

  
  char param_string[40]; // "file=" (5) + filename (up to 31) + null terminator
  snprintf(param_string, sizeof(param_string), "file=%s", agama_potential_filename);

  /* Create the AGAMA potential */
  potential->agama_potential = agama_createPotential(param_string);


  /* Read the other parameters of the model */
  potential->timestep_mult = parser_get_opt_param_float(
      parameter_file, "AGAMAPotential:timestep_mult", FLT_MAX);
}

/**
 * @brief Prints the properties of the external potential to stdout.
 *
 * @param  potential The external potential properties.
 */
static INLINE void potential_print_backend(
    const struct external_potential *potential) {

  message(
      "External potential is 'AGAMA' with properties (x,y,z) = (%e, %e, "
      "%e), timestep multiplier = %e.",
      potential->x[0], potential->x[1], potential->x[2],
      potential->timestep_mult);
}

#endif /* SWIFT_POTENTIAL_AGAMA_H */
