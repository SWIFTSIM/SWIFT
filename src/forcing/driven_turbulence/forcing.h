/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Yolan Uyttenhove (yolan.uyttenhove@ugent.be)
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
#ifndef SWIFT_FORCING_DRIVEN_TURBULENCE_H
#define SWIFT_FORCING_DRIVEN_TURBULENCE_H

/* Config parameters. */
#include <config.h>

/* Standard includes. */
#include <float.h>
#include <gsl.h>
#include <stddef.h>

/* Local includes. */
#include "error.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"

/**
 * @brief Forcing Term Properties
 */
struct forcing_terms {
  // For some small wavenumbers (large scales) store f_0 (amplitudes) for each
  // wavevector. These must be adjusted to yield the desired mach number.
  // Default: Paraboloid in fourier space, |k| in (1, 3)
  size_t wave_vectors_count;
  float *wave_vectors;

  // Store the current f (vector amplitude) for each wavevector
  float *forcing_amplitudes;

  float *f_0;

  // Store forcing parameter eta to select mixture of purely solenoidal or
  // compressive forcing power (default 0.5)
  float zeta;

  // Store inverse characteristic timescale/turnover time T
  float inv_forcing_timescale;

  // Store seeded rng
  gsl_rng *rng;
};

/**
 * @brief Computes the forcing terms.
 *
 * We do nothing in this 'none' scheme.
 *
 * @param time The current time.
 * @param terms The properties of the forcing terms.
 * @param s The #space we act on.
 * @param phys_const The physical constants in internal units.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void forcing_terms_apply(
    const double time, const struct forcing_terms* terms, const struct space* s,
    const struct phys_const* phys_const, struct part* p, struct xpart* xp) {
  /* Either update velocity directly OR update a_hydro */
  /* NOTE: for moving mesh, a_hydro is currently not used/existing */

  /* Compute fourier series with coefficients f at location of particle */
}

/**
 * @brief Updates the forcing terms.
 *
 * We do nothing in this 'none' scheme.
 *
 * @param time The current time.
 * @param timestep The timestep to the current time.
 * @param terms The properties of the forcing terms.
 * @param s The #space we act on.
 * @param us The current internal system of units
 * @param phys_const The physical constants in internal units.
 */
__attribute__((always_inline)) INLINE static void forcing_terms_update(
    const double time, const double timestep, const struct forcing_terms* terms,
    const struct space* s, const struct unit_system* us,
    const struct phys_const* phys_const) {
  // Ornstein-Uhlenbeck (OU) process. For each wavevector:
  // - subtract f * dt / T of f (exponentially decaying autocorrelation)
  // - add random "diffusion" term: f_0 * projection(eta) \dot W(dt)
  //   with W(dt) a gaussian random increment to the vector field (Wiener
  //   process).
  for (int i = 0; i < terms->wave_vectors_count; i++) {
    float *k = &terms->wave_vectors[3 * i];
    float *f = &terms->forcing_amplitudes[3 * i];
    float f_0 = terms->f_0[i];

    float df[3] = {0.f, 0.f, 0.f};

    // Exponentially decaying autocorrelation
    float diffusion_fac = timestep * terms->inv_forcing_timescale;
    df[0] -= f[0] * diffusion_fac;
    df[1] -= f[1] * diffusion_fac;
    df[2] -= f[2] * diffusion_fac;

    // random diffusion term TODO: needs normalizing/rescaling?
    float dW[3];
    dW[0] = gsl_ran_gaussian(terms->rng, timestep);
    dW[1] = gsl_ran_gaussian(terms->rng, timestep);
    dW[2] = gsl_ran_gaussian(terms->rng, timestep);
    // Project the random increment and add to df
    df[0] += f_0 * terms->zeta * dW[0];
    df[1] += f_0 * terms->zeta * dW[1];
    df[2] += f_0 * terms->zeta * dW[2];

    float k_sq = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
    float fac = f_0 * (1.f - 2.f * terms->zeta) / k_sq;
    for (int ii = 0; ii < 3; ii++) {
      df[0] += fac * k[0] * k[ii];
      df[1] += fac * k[1] * k[ii];
      df[2] += fac * k[2] * k[ii];
    }

    // Update the forcing amplitudes
    f[0] += df[0];
    f[1] += df[1];
    f[2] += df[2];
  }
}

/**
 * @brief Computes the time-step condition due to the forcing terms.
 *
 * Nothing to do here. --> Return FLT_MAX.
 *
 * @param time The current time.
 * @param terms The properties of the forcing terms.
 * @param phys_const The physical constants in internal units.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static float forcing_terms_timestep(
    double time, const struct forcing_terms* terms,
    const struct phys_const* phys_const, const struct part* p,
    const struct xpart* xp) {

  // Maybe limit based on T?

  /* No time-step size limit */
  return FLT_MAX;
}

/**
 * @brief Prints the properties of the forcing terms to stdout.
 *
 * @param terms The #forcing_terms properties of the run.
 */
static INLINE void forcing_terms_print(const struct forcing_terms* terms) {
  // TODO
  message("Forcing terms is 'No forcing terms'.");
}

/**
 * @brief Initialises the forcing term properties
 *
 * Nothing to do here.
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param s The #space object.
 * @param terms The forcing term properties to initialize
 */
static INLINE void forcing_terms_init(struct swift_params* parameter_file,
                                      const struct phys_const* phys_const,
                                      const struct unit_system* us,
                                      const struct space* s,
                                      struct forcing_terms* terms) {

  /* Nothing to do here */
}

#endif /* SWIFT_FORCING_DRIVEN_TURBULENCE_H */
