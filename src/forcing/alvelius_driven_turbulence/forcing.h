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
#include <gsl/gsl_rng.h>
#include <stddef.h>
#include <strings.h>

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
  float* wave_vectors;

  // Store the unit vectors that are mutually orthogonal and are both
  // orthogonal to the corresponding wave_vector
  float* e1;
  float* e2;

  // Store the current f (vector amplitude) for each wavevector
  float* amplitudes_real;
  float* amplitudes_imaginary;

  // The forcing for each mode
  float* kforce;

  /*! @brief Driving time step (in s). */
  double _time_step;

  /*! @brief Number of driving steps since the start of the simulation. */
  uint_fast32_t _number_of_driving_steps;

  // Store seeded rng
  gsl_rng* rng;
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
  double force[3] = {0., 0., 0.};
  for (unsigned int i = 0; i < terms->wave_vectors_count; i++) {
    const double angle = 2. * M_PI *
                         (terms->wave_vectors[3 * i] * p->x[0] +
                          terms->wave_vectors[3 * i + 1] * p->x[1] +
                          terms->wave_vectors[3 * i + 2] * p->x[2]);
    const double cos_angle = cos(angle);
    const double sin_angle = sin(angle);

    const float* fr = &terms->amplitudes_real[3 * i];
    const float* fi = &terms->amplitudes_imaginary[3 * i];
    force[0] += fr[0] * cos_angle - fi[0] * sin_angle;
    force[1] += fr[1] * cos_angle - fi[1] * sin_angle;
    force[2] += fr[2] * cos_angle - fi[2] * sin_angle;
  }
  p->a_hydro[0] += force[0];
  p->a_hydro[1] += force[1];
  p->a_hydro[2] += force[2];
}

/**
 * @brief Function gets the real and imaginary parts of the amplitudes
 * Aran and Bran of the unit vector e1 and e2, respectively, as in Eq. 11.
 * Alvelius (1999).
 *
 * @param RandGen Random number generator.
 * @param RealRand Output array to store the real parts of the amplitudes of
 * the forcing.
 * @param ImRand Output array to store the imaginary parts of the amplitudes
 * of the forcing.
 */
static void get_random_factors(gsl_rng* rng, double* RealRand, double* ImRand) {

  const double phi = 2. * M_PI * gsl_rng_uniform(rng);
  const double ga = sin(phi);
  const double gb = cos(phi);
  const double theta1 = 2. * M_PI * gsl_rng_uniform(rng);
  const double theta2 = 2. * M_PI * gsl_rng_uniform(rng);
  RealRand[0] = cos(theta1) * ga;
  ImRand[0] = sin(theta1) * ga;
  RealRand[1] = cos(theta2) * gb;
  ImRand[1] = sin(theta2) * gb;
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
    const double time, const double timestep, struct forcing_terms* terms,
    const struct space* s, const struct unit_system* us,
    const struct phys_const* phys_const) {
  /* Reset amplitudes */
  bzero(terms->amplitudes_real,
        3 * terms->wave_vectors_count * sizeof(*terms->amplitudes_real));
  bzero(terms->amplitudes_imaginary,
        3 * terms->wave_vectors_count * sizeof(*terms->amplitudes_imaginary));

  /* Compute combined contribution of the required number of driving steps */
  while (terms->_number_of_driving_steps * terms->_time_step < time) {
    for (unsigned int i = 0; i < terms->wave_vectors_count; ++i) {
      double RealRand[2];
      double ImRand[2];
      get_random_factors(terms->rng, RealRand, ImRand);

      for (int j = 0; j < 3; j++) {
        terms->amplitudes_real[i] +=
            terms->kforce[i] * terms->e1[3 * i + j] * RealRand[0] +
            terms->kforce[i] * terms->e2[3 * i + j] * RealRand[1];
        terms->amplitudes_imaginary[i] +=
            terms->kforce[i] * terms->e1[3 * i + j] * ImRand[0] +
            terms->kforce[i] * terms->e2[3 * i + j] * ImRand[1];
      }
    }
    terms->_number_of_driving_steps++;
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

  /* No time-step size limit */
  return FLT_MAX;
}

/**
 * @brief Prints the properties of the forcing terms to stdout.
 *
 * @param terms The #forcing_terms properties of the run.
 */
static INLINE void forcing_terms_print(const struct forcing_terms* terms) {
  message("Forcing terms is 'Alvelius driven turbulence'.");
  message("Forcing 'Alvelius driven turbulence' with %lu modes",
          terms->wave_vectors_count);
}

/**
 * @brief Initialises the forcing term properties
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
  /* Read parameters */
  const double k_min = parser_get_opt_param_double(
      parameter_file, "AlveliusTurbulenceForcing:minimum_wave_number", 1.);
  const double k_max = parser_get_opt_param_double(
      parameter_file, "AlveliusTurbulenceForcing:maximum_wave_number", 3.);
  const double k_peak = parser_get_opt_param_double(
      parameter_file, "AlveliusTurbulenceForcing:peak_forcing_wave_number",
      2.5);
  const double concentration_factor = parser_get_opt_param_double(
      parameter_file, "AlveliusTurbulenceForcing:concentration_factor", 0.2);
  const double forcing_power_unit_factor =
      us->UnitTime_in_cgs * us->UnitTime_in_cgs * us->UnitTime_in_cgs /
      (us->UnitLength_in_cgs * us->UnitLength_in_cgs);
  const double forcing_power =
      forcing_power_unit_factor *
      parser_get_opt_param_double(parameter_file,
                                  "AlveliusTurbulenceForcing:forcing_power_cgs",
                                  2.717e-2 /* cm^2 s^-3 */);
  double time_step =
      parser_get_opt_param_double(parameter_file,
                                  "AlveliusTurbulenceForcing:time_step_cgs",
                                  1.519e6 /* s */) /
      us->UnitTime_in_cgs;
  const int seed = parser_get_opt_param_int(
      parameter_file, "AlveliusTurbulenceForcing:random_seed", 42);

  /* Get boxsize */
  const double Linv = 1. / max3(s->dim[0], s->dim[1], s->dim[2]);

  /* The force spectrum here prescribed is  Gaussian in shape:
   * F(k) = amplitude*exp^((k-k_peak)^2/concentration_factor)
   *        amplitude*gaussian_exponential
   */
  double spectra_sum = 0.;
  const double cinv = 1. / (concentration_factor * concentration_factor);

  /* Construct the wavevectors and the two orthogonal directions along which the
   * force will be applied */
  size_t number_of_modes = 0;
  size_t wave_vectors_size = 10;
  float* wave_vectors = swift_malloc(
      "Forcing terms", 3 * wave_vectors_size * sizeof(*wave_vectors));
  float* e1_arr =
      swift_malloc("Forcing terms", 3 * wave_vectors_size * sizeof(*e1_arr));
  float* e2_arr =
      swift_malloc("Forcing terms", 3 * wave_vectors_size * sizeof(*e2_arr));
  float* kforce =
      swift_malloc("Forcing terms", wave_vectors_size * sizeof(*kforce));
  for (double k1 = 0.; k1 <= k_max; k1 += 1.) {
    double k2start;
    if (k1 == 0.) {
      k2start = 0.;
    } else {
      k2start = -k_max;
    }
    for (double k2 = k2start; k2 <= k_max; k2 += 1.) {
      double k3start;
      if (k1 == 0. && k2 == 0.) {
        k3start = 0.;
      } else {
        k3start = -k_max;
      }
      for (double k3 = k3start; k3 <= k_max; k3 += 1.) {
        const double pwrk1 = k1 * k1;
        const double pwrk2 = k2 * k2;
        const double pwrk3 = k3 * k3;
        const double kk = pwrk1 + pwrk2 + pwrk3;
        const double k = sqrt(kk);
        if (k <= k_max && k >= k_min) {
          const double kdiff = (k - k_peak);
          const double sqrtk12 = sqrt(pwrk1 + pwrk2);
          const double invkk = 1. / kk;
          const double invk = 1. / k;

          double e1[3], e2[3];
          if (sqrtk12 > 0.) {
            const double invsqrtk12 = 1. / sqrtk12;
            e1[0] = k2 * invsqrtk12;
            e1[1] = -k1 * invsqrtk12;
            e1[2] = 0.;
            e2[0] = k1 * k3 * invsqrtk12 * invk;
            e2[1] = k2 * k3 * invsqrtk12 * invk;
            e2[2] = -sqrtk12 * invk;
          } else {
            const double sqrtk13 = sqrt(pwrk1 + pwrk3);
            const double invsqrtk13 = 1. / sqrtk13;
            e1[0] = -k3 * invsqrtk13;
            e1[1] = 0.;
            e1[2] = k1 * invsqrtk13;
            e2[0] = k1 * k2 * invsqrtk13 * invk;
            e2[1] = -sqrtk13 * invk;
            e2[2] = k2 * k3 * invsqrtk13 * invk;
          }
          assert(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2] <= 1.1);
          assert(e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2] <= 1.1);

          const double gaussian_spectra = exp(-kdiff * kdiff * cinv) * invkk;
          spectra_sum += gaussian_spectra;

          /* Do we still have space in the arrays? */
          if (number_of_modes == wave_vectors_size) {
            wave_vectors_size <<= 1;
            wave_vectors =
                swift_realloc("Forcing terms", wave_vectors,
                              3 * wave_vectors_size * sizeof(*wave_vectors));
            e1_arr = swift_realloc("Forcing terms", e1_arr,
                                   3 * wave_vectors_size * sizeof(*e1_arr));
            e2_arr = swift_realloc("Forcing terms", e2_arr,
                                   3 * wave_vectors_size * sizeof(*e2_arr));
            kforce = swift_realloc("Forcing terms", kforce,
                                   wave_vectors_size * sizeof(*kforce));
          }

          /* Append to arrays */
          wave_vectors[3 * number_of_modes] = k1 * Linv;
          wave_vectors[3 * number_of_modes + 1] = k2 * Linv;
          wave_vectors[3 * number_of_modes + 2] = k3 * Linv;
          e1_arr[3 * number_of_modes] = e1[0];
          e1_arr[3 * number_of_modes + 1] = e1[1];
          e1_arr[3 * number_of_modes + 2] = e1[2];
          e2_arr[3 * number_of_modes] = e2[0];
          e2_arr[3 * number_of_modes + 1] = e2[1];
          e2_arr[3 * number_of_modes + 2] = e2[2];
          kforce[number_of_modes] = gaussian_spectra;
          number_of_modes++;
        }
      }
    }
  }

  /* Obtain full expression for the forcing amplitude */
  const float norm = forcing_power / (spectra_sum * time_step);
  for (size_t i = 0; i < number_of_modes; ++i) {
    kforce[i] *= norm;
    kforce[i] = sqrtf(kforce[i]);
  }

  /* Finally initialize the struct */
  terms->wave_vectors_count = number_of_modes;
  terms->wave_vectors = wave_vectors;
  terms->e1 = e1_arr;
  terms->e2 = e2_arr;
  terms->kforce = kforce;
  terms->amplitudes_real = swift_calloc("Forcing terms", 3 * number_of_modes,
                                        sizeof *terms->amplitudes_real);
  terms->amplitudes_imaginary =
      swift_calloc("Forcing terms", 3 * number_of_modes,
                   sizeof *terms->amplitudes_imaginary);
  terms->_time_step = time_step;
  terms->_number_of_driving_steps = 0;
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, seed);
  terms->rng = rng;
}

/**
 * @brief Frees any memory allocated for the forcing terms
 *
 * @param terms The forcing term properties to clean
 */
static INLINE void forcing_terms_clean(struct forcing_terms* terms) {
  swift_free("Forcing terms", terms->wave_vectors);
  swift_free("Forcing terms", terms->e1);
  swift_free("Forcing terms", terms->e2);
  swift_free("Forcing terms", terms->kforce);
  swift_free("Forcing terms", terms->amplitudes_real);
  swift_free("Forcing terms", terms->amplitudes_imaginary);
  gsl_rng_free(terms->rng);
}

#endif /* SWIFT_FORCING_DRIVEN_TURBULENCE_H */
