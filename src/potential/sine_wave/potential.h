/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_SINE_WAVE_H
#define SWIFT_SINE_WAVE_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <math.h>

/* Local includes. */
#include "const.h"
#include "error.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"

/**
 * @brief External Potential Properties - Sine wave case
 */
struct external_potential {

  /*! Amplitude of the sine wave. */
  double amplitude;

  /*! Growth time of the potential. */
  double growth_time;

  /*! Time-step limiting factor. */
  double timestep_limit;
};

/**
 * @brief Computes the time-step from the acceleration due to a sine wave.
 *
 * @param time The current time.
 * @param potential The properties of the potential.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static float external_gravity_timestep(
    double time, const struct external_potential* restrict potential,
    const struct phys_const* restrict phys_const,
    const struct gpart* restrict g) {

  return potential->timestep_limit;
}

/**
 * @brief Computes the gravitational acceleration along x given by the sine
 * wave.
 *
 * @param time The current time in internal units.
 * @param potential The properties of the potential.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static void external_gravity_acceleration(
    double time, const struct external_potential* restrict potential,
    const struct phys_const* restrict phys_const, struct gpart* restrict g) {

  float Acorr = 1.;

  if (time < potential->growth_time) {
    Acorr = time / potential->growth_time;
  }

  g->a_grav[0] = potential->amplitude * Acorr * sin(2. * M_PI * g->x[0]) /
                 phys_const->const_newton_G;
}

/**
 * @brief Computes the gravitational potential energy of a particle in the
 * sine wave.
 *
 * @param time The current time.
 * @param potential The #external_potential used in the run.
 * @param phys_const Physical constants in internal units.
 * @param gp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
external_gravity_get_potential_energy(
    double time, const struct external_potential* potential,
    const struct phys_const* const phys_const, const struct gpart* gp) {

  /* this potential does not really have a potential energy */
  return 0.;
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
    const struct phys_const* phys_const, const struct unit_system* us,
    const struct space* s, struct external_potential* potential) {

  potential->amplitude =
      parser_get_param_double(parameter_file, "SineWavePotential:amplitude");
  potential->growth_time = parser_get_opt_param_double(
      parameter_file, "SineWavePotential:growth_time", 0.);
  potential->timestep_limit = parser_get_param_double(
      parameter_file, "SineWavePotential:timestep_limit");
}

/**
 * @brief Prints the properties of the external potential to stdout.
 *
 * @param  potential The external potential properties.
 */
static INLINE void potential_print_backend(
    const struct external_potential* potential) {

  message("External potential is a sine wave with amplitude %g",
          potential->amplitude);
}

#endif /* SWIFT_SINE_WAVE_H */
