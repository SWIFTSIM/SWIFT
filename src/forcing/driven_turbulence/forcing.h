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
#ifndef SWIFT_FORCING_NONE_H
#define SWIFT_FORCING_NONE_H

/* Config parameters. */
#include <config.h>

/* Standard includes. */
#include <float.h>

/* Library includes */
#include <gsl/gsl_rng.h>

/* Local includes. */
#include "dimension.h"
#include "error.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"

enum amplitude_spectrum {
  constant,
  parabolic, /*! Federrath et al. 2010 */
  kolmogorov,
  quadratic,
};

/**
 * @brief Forcing Term Properties
 */
struct forcing_terms {

  /*! Shape of the power spectrum */
  enum amplitude_spectrum shape;

  /*! Number of modes considered */
  int num_modes;

  /*! Minimal mode */
  double k_min;

  /*! Maximal mode */
  double k_max;

  double time_frequency;  // ST_DtFreq

  double decay_time;  // ST_Decay

  /*! Ornstein-Uhlenbeck varience */
  double variance;  // StOUVar

  double solenoid_weight;  // ST_SolWeight

  double solenoid_weight_norm;  // StSolWeightNorm

  double amplitude_factor;  // ST_AmplFac

  /*! Modes considered */
  double *modes;

  /*! Amplitudes of the modes */
  double *amplitudes;

  /*! Ornstein-Uhlenbeck phases */
  double *phases;

  /*! Complex amplitudes */
  double *A_k_a;
  double *A_k_b;

  /*! Random number generator */
  gsl_rng *rng;
};

/**
 * @brief Computes the hydrodynamic forcing terms.
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
__attribute__((always_inline)) INLINE static void forcing_hydro_terms_apply(
    const double time, const struct forcing_terms *terms, const struct space *s,
    const struct phys_const *phys_const, struct part *p, struct xpart *xp) {

  double fx = 0.;
  double fy = 0.;
  double fz = 0.;
  const double rho = hydro_get_comoving_density(p);

  for (int m = 0; m < terms->num_modes; ++m) {

    const double kxx = terms->modes[3 * m + 0] * p->x[0];
    const double kyy = terms->modes[3 * m + 1] * p->x[1];
    const double kzz = terms->modes[3 * m + 2] * p->x[2];
    const double k_dot_x = kxx + kyy + kzz;

    const double A = terms->amplitudes[m];

    double real, imag;
    sincos(k_dot_x, &real, &imag);

    fx += A * (terms->A_k_a[3 * m + 0] * real - terms->A_k_b[3 * m + 0] * imag);
    fy += A * (terms->A_k_a[3 * m + 1] * real - terms->A_k_b[3 * m + 1] * imag);
    fz += A * (terms->A_k_a[3 * m + 2] * real - terms->A_k_b[3 * m + 2] * imag);
  }

  /* Turbulent force */
  fx *= 2. * terms->amplitude_factor * terms->solenoid_weight_norm;
  fy *= 2. * terms->amplitude_factor * terms->solenoid_weight_norm;
  fz *= 2. * terms->amplitude_factor * terms->solenoid_weight_norm;

  p->a_hydro[0] += fx / rho;
  p->a_hydro[1] += fy / rho;
  p->a_hydro[2] += fz / rho;
}

/**
 * @brief Computes the gravitational forcing terms.
 *
 * We do nothing in this 'none' scheme.
 *
 * @param id The particle ID.
 * @param terms The properties of the forcing terms.
 * @param gp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static void forcing_grav_terms_apply(
    const long long id, const struct forcing_terms *terms, struct gpart *gp) {
  /* Nothing to do here */
}

/**
 * @brief Sets the forcing of gparts prior to drift.
 *
 * @param id The particle ID.
 * @param terms The properties of the forcing terms.
 * @param gp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static void forcing_gpart_drift_apply(
    const long long id, const struct forcing_terms *terms, struct gpart *gp) {}

/**
 * @brief Sets the forcing of parts prior to drift.
 *
 * @param id The particle ID.
 * @param terms The properties of the forcing terms.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void forcing_part_drift_apply(
    const long long id, const struct forcing_terms *terms, struct part *p,
    struct xpart *xp) {}

/**
 * @brief Sets the forcing of sparts prior to drift.
 *
 * @param id The particle ID.
 * @param terms The properties of the forcing terms.
 * @param sp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static void forcing_spart_drift_apply(
    const long long id, const struct forcing_terms *terms, struct spart *sp) {}

/**
 * @brief Sets the forcing of bparts prior to drift.
 *
 * @param id The particle ID.
 * @param terms The properties of the forcing terms.
 * @param bp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static void forcing_bpart_drift_apply(
    const long long id, const struct forcing_terms *terms, struct bpart *bp) {}

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
    double time, const struct forcing_terms *terms, const struct space *s,
    const struct phys_const *phys_const, const struct part *p,
    const struct xpart *xp) {

  /* No time-step size limit */
  return FLT_MAX;
}

static INLINE void forcing_terms_init_random_sequence(
    struct forcing_terms *terms) {

  for (int i = 0; i < 6 * terms->num_modes; i++) {

    /* Random normally-distributed number */
    const double r0 = gsl_rng_uniform(terms->rng);
    const double r1 = gsl_rng_uniform(terms->rng);
    const double n = sqrt(2. * log(1. / r0)) * cos(2. * M_PI * r1);

    terms->phases[i] = n * terms->variance;
  }
}

static INLINE void forcing_terms_update_random_sequence(
    struct forcing_terms *terms) {

  const double damping = exp(-terms->time_frequency / terms->decay_time);

  const double term2 = terms->variance * sqrt(1. - damping * damping);

  for (int i = 0; i < 6 * terms->num_modes; i++) {

    /* Random normally-distributed number */
    const double r0 = gsl_rng_uniform(terms->rng);
    const double r1 = gsl_rng_uniform(terms->rng);
    const double n = sqrt(2. * log(1. / r0)) * cos(2. * M_PI * r1);

    terms->phases[i] = terms->phases[i] * damping + term2 * n;
  }
}

static INLINE void forcing_terms_calculate_phases(struct forcing_terms *terms) {

  const double solenoid_weight = terms->solenoid_weight;

  for (int i = 0; i < terms->num_modes; ++i) {

    double k_a = 0., k_b = 0., k_k = 0.;

    for (int j = 0; j < hydro_dimension_integer; j++) {

      k_k += terms->modes[3 * i + j] * terms->modes[3 * i + j];
      k_a += terms->modes[3 * i + j] * terms->phases[6 * i + 2 * j + 1];
      k_b += terms->modes[3 * i + j] * terms->phases[6 * i + 2 * j + 0];
    }

    for (int j = 0; j < hydro_dimension_integer; j++) {

      const double div_a = terms->modes[3 * i + j] * k_a / k_k;
      const double div_b = terms->modes[3 * i + j] * k_b / k_k;
      const double curl_a = terms->phases[6 * i + 2 * j + 0] - div_b;
      const double curl_b = terms->phases[6 * i + 2 * j + 0] - div_a;

      terms->A_k_a[3 * i + j] =
          solenoid_weight * curl_a + (1. - solenoid_weight) * div_b;
      terms->A_k_b[3 * i + j] =
          solenoid_weight * curl_b + (1. - solenoid_weight) * div_a;
    }
  }
}

/**
 * @brief updates the forcing terms
 *
 * Nothing to do here
 *
 * @param terms The #forcing_terms properties of the run
 * @param time_old The previous system time
 */
INLINE static void forcing_update(struct forcing_terms *terms,
                                  const double time_old) {

  // CHECK IF TIME HAS ELAPSED

  forcing_terms_init_random_sequence(terms);
  forcing_terms_calculate_phases(terms);
}

/**
 * @brief Prints the properties of the forcing terms to stdout.
 *
 * @param terms The #forcing_terms properties of the run.
 */
static INLINE void forcing_terms_print(const struct forcing_terms *terms) {

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
static INLINE void forcing_terms_init(struct swift_params *parameter_file,
                                      const struct phys_const *phys_const,
                                      const struct unit_system *us,
                                      const struct space *s,
                                      struct forcing_terms *terms) {

  /* Read in parameter files ---------------------------------*/

  int seed = 0;
  double energy = 1;
  double decay = 1;
  double weight = 1;  // ST_SolWeight

  terms->amplitude_factor = 1.0;  // ST_AmplFac
  terms->variance = sqrt(energy / decay);
  terms->solenoid_weight = weight;
  terms->solenoid_weight_norm =  // StSolWeightNorm
      sqrt(3.0 / hydro_dimension) * sqrt(3.0) * 1.0 /
      sqrt(1.0 - 2.0 * weight + hydro_dimension * weight * weight);

  /* Prepare the random number generator ---------------------*/

  /* Note: we use GSL and not our mechanism to match the choices
   * of Pakmor+2026 */
  terms->rng = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(terms->rng, seed);

  /* Compute number of modes -------------------------------- */

  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};

  const double a_min = 0.;
  const double k_min = terms->k_min;
  const double k_max = terms->k_max;
  const double k_mid = 0.5 * (k_max + k_min);
  const enum amplitude_spectrum shape = terms->shape;

#if HYDRO_DIMENSION_1D
  const int ikx_max = dim[0] * k_max / 2. / M_PI;
  const int iky_max = 0;
  const int ikz_max = 0;
#elif HYDRO_DIMENSION_2D
  const int ikx_max = dim[0] * k_max / 2. / M_PI;
  const int iky_max = dim[1] * k_max / 2. / M_PI;
  const int ikz_max = 0;
#else
  const int ikx_max = dim[0] * k_max / 2. / M_PI;
  const int iky_max = dim[1] * k_max / 2. / M_PI;
  const int ikz_max = dim[2] * k_max / 2. / M_PI;
#endif

  int num_modes = 0;
  for (int ikx = 0; ikx <= ikx_max; ikx++) {
    for (int iky = 0; iky <= iky_max; iky++) {
      for (int ikz = 0; ikz <= ikz_max; ikz++) {

        const double kx = 2. * M_PI * ikx / dim[0];
        const double ky = 2. * M_PI * iky / dim[1];
        const double kz = 2. * M_PI * ikz / dim[2];
        const double k = sqrt(kx * kx + ky * ky + kz * kz);

        if (k >= k_min && k <= k_max) {
#if HYDRO_DIMENSION_1D
          num_modes += 1;
#elif HYDRO_DIMENSION_2D
          num_modes += 2;
#else
          num_modes += 4;
#endif
        }
      }
    }
  }

  /* Allocate modes ----------------------------------------- */

  terms->modes =
      (double *)swift_malloc("forcing_modes", num_modes * 3 * sizeof(double));
  if (terms->modes == NULL) error("Error allocating forcing modes array");

  terms->amplitudes =
      (double *)swift_malloc("forcing_amplitudes", num_modes * sizeof(double));
  if (terms->amplitudes == NULL)
    error("Error allocating forcing amplitudes array");

  terms->phases =
      (double *)swift_malloc("forcing_phases", num_modes * 6 * sizeof(double));
  if (terms->phases == NULL) error("Error allocating forcing phases array");

  terms->A_k_a =
      (double *)swift_malloc("forcing_Aka", num_modes * 3 * sizeof(double));
  if (terms->A_k_a == NULL)
    error("Error allocating forcing complex amplitude array");

  terms->A_k_b =
      (double *)swift_malloc("forcing_Akb", num_modes * 3 * sizeof(double));
  if (terms->A_k_b == NULL)
    error("Error allocating forcing complex amplitude array");

  /* Compute modes ------------------------------------------ */

  num_modes = 0;
  for (int ikx = 0; ikx <= ikx_max; ikx++) {
    for (int iky = 0; iky <= iky_max; iky++) {
      for (int ikz = 0; ikz <= ikz_max; ikz++) {

        const double kx = 2. * M_PI * ikx / dim[0];
        const double ky = 2. * M_PI * iky / dim[1];
        const double kz = 2. * M_PI * ikz / dim[2];
        const double k = sqrt(kx * kx + ky * ky + kz * kz);

        if (k >= k_min && k <= k_max) {

          double amplitude = 0.;
          switch (shape) {
            case constant:
              amplitude = 1.;
              break;
            case parabolic:
              amplitude = 4.0 * (a_min - 1.0) /
                              ((k_max - k_min) * (k_max - k_min)) *
                              ((k - k_mid) * (k - k_mid)) +
                          1.0;
              break;
            case kolmogorov:
              amplitude = pow(k_min / k, 5. / 3.);
              break;
            case quadratic:
              amplitude = pow(k_min / k, 2.);
              break;
          }

          terms->amplitudes[num_modes] = amplitude;
          terms->modes[3 * num_modes + 0] = kx;
          terms->modes[3 * num_modes + 1] = ky;
          terms->modes[3 * num_modes + 2] = kz;
          num_modes++;
#if defined(HYDRO_DIMENSION_2D) || defined(HYDRO_DIMENSION_3D)
          terms->amplitudes[num_modes] = amplitude;
          terms->modes[3 * num_modes + 0] = kx;
          terms->modes[3 * num_modes + 1] = -ky;
          terms->modes[3 * num_modes + 2] = kz;
          num_modes++;
#endif
#if defined(HYDRO_DIMENSION_3D)
          terms->amplitudes[num_modes] = amplitude;
          terms->modes[3 * num_modes + 0] = kx;
          terms->modes[3 * num_modes + 1] = ky;
          terms->modes[3 * num_modes + 2] = -kz;
          num_modes++;

          terms->amplitudes[num_modes] = amplitude;
          terms->modes[3 * num_modes + 0] = kx;
          terms->modes[3 * num_modes + 1] = -ky;
          terms->modes[3 * num_modes + 2] = -kz;
          num_modes++;
#endif
        }
      }
    }
  }

  terms->num_modes = num_modes;

  /* Initialise the Ornstein-Uhlenbeck phases ------------------ */

  forcing_terms_init_random_sequence(terms);

  /* Initialise the complex phases ------------------------------ */

  forcing_terms_calculate_phases(terms);
}

/**
 * @brief Clean-up the memory allocated for the forcing routine
 *
 * Nothing to do here
 *
 * @param terms The forcing term properties
 */
static INLINE void forcing_terms_clean(struct forcing_terms *terms) {

  swift_free("forcing_modes", terms->modes);
  swift_free("forcing_amplitudes", terms->amplitudes);
  swift_free("forcing_phases", terms->phases);
  swift_free("forcing_Aka", terms->A_k_a);
  swift_free("forcing_Akb", terms->A_k_b);
}

#endif /* SWIFT_FORCING_NONE_H */
