/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *               2021 Nina Sartorio (sartorio.nina@gmail.com)
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
#ifndef SWIFT_TURBULENCE_H
#define SWIFT_TURBULENCE_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include "engine.h"

#include <float.h>
#include <gsl/gsl_rng.h>

#ifdef TURBULENCE_DRIVING_ALVELIUS

/**
 * @brief Properties of the Alvelius turbulence driving.
 */
struct turbulence_driving {

  /* Random number generator. */
  gsl_rng* random_generator;

  /* Driving time step. */
  double dt;

  /* Number of driving steps already taken. */
  int number_of_steps;

  /* Number of modes in Fourier space. */
  int number_of_modes;

  /* Wave vectors in Fourier space. */
  double* k;

  /* Amplitudes of the forcing. */
  double* amplitudes;

  /* Unit vectors of the forcing. */
  double* unit_vectors;

  /* Forcing of each mode. */
  double* forcing;
};

#else /* TURBULENCE_DRIVING_NONE */

/**
 * @brief Empty turbulence driving structure.
 */
struct turbulence_driving {};

#endif /* TURBULENCE_DRIVING */

void turbulence_accelerate(struct part* restrict p, struct xpart* restrict xp,
                           const struct turbulence_driving* restrict
                               turbulence);
void turbulence_init_backend(struct swift_params* parameter_file,
                             const struct phys_const* phys_const,
                             const struct unit_system* us,
                             const struct space* s,
                             struct turbulence_driving* turbulence);
void turbulence_print_backend(const struct turbulence_driving* turbulence);
void turbulence_update(struct engine* restrict e);

#endif /* SWIFT_TURBULENCE_H */
