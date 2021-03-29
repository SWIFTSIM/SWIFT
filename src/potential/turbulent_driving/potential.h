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
#ifndef SWIFT_POTENTIAL_TURBULENT_DRIVING_H
#define SWIFT_POTENTIAL_TURBULENT_DRIVING_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <gsl/gsl_rng.h>

/* Local includes. */
#include "error.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"

/**
 * @brief External Potential Properties
 */
struct external_potential {

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

/**
 * @brief Computes the time-step due to the acceleration.
 *
 * @param time The current time.
 * @param potential The properties of the externa potential.
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
 * @brief Computes the gravitational acceleration due to nothing
 *
 * We do nothing.
 *
 * @param time The current time.
 * @param potential The proerties of the external potential.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static void external_gravity_acceleration(
    double time, const struct external_potential* restrict potential,
    const struct phys_const* restrict phys_const, struct gpart* restrict g) {

  const double* x = g->x;

  double force[3] = {0., 0., 0.};
  for (int ik = 0; ik < potential->number_of_modes; ++ik) {
    const double fr[3] = {potential->amplitudes[6 * ik + 0],
                          potential->amplitudes[6 * ik + 1],
                          potential->amplitudes[6 * ik + 2]};
    const double fi[3] = {potential->amplitudes[6 * ik + 3],
                          potential->amplitudes[6 * ik + 4],
                          potential->amplitudes[6 * ik + 5]};
    const double k[3] = {potential->k[3 * ik + 0], potential->k[3 * ik + 1],
                         potential->k[3 * ik + 2]};

    const double cosx = cos(2. * M_PI * k[0] * x[0]);
    const double cosy = cos(2. * M_PI * k[1] * x[1]);
    const double cosz = cos(2. * M_PI * k[2] * x[2]);
    const double sinx = sin(2. * M_PI * k[0] * x[0]);
    const double siny = sin(2. * M_PI * k[1] * x[1]);
    const double sinz = sin(2. * M_PI * k[2] * x[2]);

    const double cosyz = cosy * cosz - siny * sinz;
    const double sinyz = siny * cosz + cosy * sinz;

    const double cosxyz = cosx * cosyz - sinx * sinyz;
    const double sinxyz = sinx * cosyz + cosx * sinyz;

    force[0] += fr[0] * cosxyz - fi[0] * sinxyz;
    force[1] += fr[1] * cosxyz - fi[1] * sinxyz;
    force[2] += fr[2] * cosxyz - fi[2] * sinxyz;
  }

  g->a_grav[0] = force[0];
  g->a_grav[1] = force[1];
  g->a_grav[2] = force[2];
}

/**
 * @brief Computes the gravitational potential energy due to nothing.
 *
 * We return 0.
 *
 * @param time The current time.
 * @param potential The #external_potential used in the run.
 * @param phys_const Physical constants in internal units.
 * @param g Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
external_gravity_get_potential_energy(
    double time, const struct external_potential* potential,
    const struct phys_const* const phys_const, const struct gpart* g) {

  return 0.f;
}

/**
 * @brief Initialises the external potential properties in the internal system
 * of units.
 *
 * Nothing to do here.
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

  /* make sure the box is a cube */
  if (s->dim[0] != s->dim[1] || s->dim[0] != s->dim[2]) {
    error("Turbulent forcing only works in a cubic box!");
  }

  /* get dimensionless parameters */
  const int random_seed = parser_get_opt_param_int(
      parameter_file, "TurbulentDrivingPotential:random_seed", 42);
  const double kmin = parser_get_opt_param_double(
      parameter_file, "TurbulentDrivingPotential:kmin", 2.);
  const double kmax = parser_get_opt_param_double(
      parameter_file, "TurbulentDrivingPotential:kmax", 3.);
  const double kforcing = parser_get_opt_param_double(
      parameter_file, "TurbulentDrivingPotential:kforcing", 2.5);
  const double concentration_factor = parser_get_opt_param_double(
      parameter_file, "TurbulentDrivingPotential:concentration_factor", 0.2);

  /* get parameters with units */
  const double power_forcing_cgs = parser_get_opt_param_double(
      parameter_file, "TurbulentDrivingPotential:power_forcing_in_cm2_per_s3",
      17.);
  const double dtfor_cgs = parser_get_opt_param_double(
      parameter_file, "TurbulentDrivingPotential:dt_forcing_in_s", 1.e6);
  const double starting_time_cgs = parser_get_opt_param_double(
      parameter_file, "TurbulentDrivingPotential:starting_time_in_s", 0.);

  /* convert units */
  const float forcing_quantity[5] = {0.0f, 2.0f, -3.0f, 0.0f, 0.0f};
  const double uf_in_cgs =
      units_general_cgs_conversion_factor(us, forcing_quantity);
  const double power_forcing = power_forcing_cgs / uf_in_cgs;
  const double ut_in_cgs = units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  const double dtfor = dtfor_cgs / ut_in_cgs;
  const double starting_time = starting_time_cgs / ut_in_cgs;

  /* pre-compute some constants */
  const double Linv = 1. / s->dim[0];
  const double cinv = 1. / (concentration_factor * concentration_factor);

  /* initialise the random number generator */
  potential->random_generator = gsl_rng_alloc(gsl_rng_ranlxd2);
  gsl_rng_set(potential->random_generator, random_seed);

  /* count the number of k-modes within the forcing shell */
  int num_modes = 0;
  for (double k1 = 0.; k1 <= kmax; k1 += 1.) {
    double k2start;
    if (k1 == 0.) {
      k2start = 0.;
    } else {
      k2start = -kmax;
    }
    for (double k2 = k2start; k2 <= kmax; k2 += 1.) {
      double k3start;
      if (k1 == 0. && k2 == 0.) {
        k3start = 0.;
      } else {
        k3start = -kmax;
      }
      for (double k3 = k3start; k3 <= kmax; k3 += 1.) {
        const double pwrk1 = k1 * k1;
        const double pwrk2 = k2 * k2;
        const double pwrk3 = k3 * k3;
        const double kk = pwrk1 + pwrk2 + pwrk3;
        const double k = sqrt(kk);
        if (k <= kmax && k >= kmin) {
          ++num_modes;
        }
      }
    }
  }

  /* allocate forcing arrays */
  potential->k =
      (double*)swift_malloc("turbulent_k", 3 * num_modes * sizeof(double));
  potential->amplitudes =
      (double*)swift_malloc("turbulent_A", 6 * num_modes * sizeof(double));
  potential->unit_vectors =
      (double*)swift_malloc("turbulent_uv", 6 * num_modes * sizeof(double));
  potential->forcing =
      (double*)swift_malloc("turbulent_f", num_modes * sizeof(double));

  /* compute the k-modes, unit vectors and forcing per k-mode */
  int kindex = 0;
  double spectrum_sum = 0.;
  for (double k1 = 0.; k1 <= kmax; k1 += 1.) {
    double k2start;
    if (k1 == 0.) {
      k2start = 0.;
    } else {
      k2start = -kmax;
    }
    for (double k2 = k2start; k2 <= kmax; k2 += 1.) {
      double k3start;
      if (k1 == 0. && k2 == 0.) {
        k3start = 0.;
      } else {
        k3start = -kmax;
      }
      for (double k3 = k3start; k3 <= kmax; k3 += 1.) {
        const double pwrk1 = k1 * k1;
        const double pwrk2 = k2 * k2;
        const double pwrk3 = k3 * k3;
        const double kk = pwrk1 + pwrk2 + pwrk3;
        const double k = sqrt(kk);
        if (k <= kmax && k >= kmin) {
          const double kdiff = (k - kforcing);
          const double sqrtk12 = sqrt(pwrk1 + pwrk2);
          const double invkk = 1. / kk;
          const double invk = 1. / k;
          if (sqrtk12 > 0.) {
            const double invsqrtk12 = 1. / sqrtk12;
            potential->unit_vectors[6 * kindex + 0] = k2 * invsqrtk12;
            potential->unit_vectors[6 * kindex + 1] = -k1 * invsqrtk12;
            potential->unit_vectors[6 * kindex + 2] = 0.;
            potential->unit_vectors[6 * kindex + 3] =
                k1 * k3 * invsqrtk12 * invk;
            potential->unit_vectors[6 * kindex + 4] =
                k2 * k3 * invsqrtk12 * invk;
            potential->unit_vectors[6 * kindex + 5] = -sqrtk12 * invk;
          } else {
            const double sqrtk13 = sqrt(pwrk1 + pwrk3);
            const double invsqrtk13 = 1. / sqrtk13;
            potential->unit_vectors[6 * kindex + 0] = -k3 * invsqrtk13;
            potential->unit_vectors[6 * kindex + 1] = 0.;
            potential->unit_vectors[6 * kindex + 2] = k1 * invsqrtk13;
            potential->unit_vectors[6 * kindex + 3] =
                k1 * k3 * invsqrtk13 * invk;
            potential->unit_vectors[6 * kindex + 4] = -sqrtk13 * invk;
            potential->unit_vectors[6 * kindex + 5] =
                k2 * k3 * invsqrtk13 * invk;
          }

          potential->k[3 * kindex + 0] = k1 * Linv;
          potential->k[3 * kindex + 1] = k2 * Linv;
          potential->k[3 * kindex + 2] = k3 * Linv;
          const double gaussian_spectrum = exp(-kdiff * kdiff * cinv) * invkk;
          spectrum_sum += gaussian_spectrum;
          potential->forcing[kindex] = gaussian_spectrum;
          ++kindex;
        }
      }
    }
  }

  /* normalise the forcing */
  const double norm = power_forcing / (spectrum_sum * dtfor);
  for (int i = 0; i < num_modes; ++i) {
    potential->forcing[i] *= norm;
    potential->forcing[i] = sqrt(potential->forcing[i]);
  }

  /* fast-forward the driving to the desired point in time */
  int num_steps = 0;
  while (num_steps * dtfor < starting_time) {
    /* 3 random numbers are generated per mode in potential_update() */
    for (int i = 0; i < 3 * num_modes; ++i) {
      gsl_rng_uniform(potential->random_generator);
    }
    ++num_steps;
  }

  /* initialise the final variables */
  potential->number_of_steps = 0;
  potential->dt = dtfor;
  potential->number_of_modes = num_modes;
}

/**
 * @brief Prints the properties of the external potential to stdout.
 *
 * @param  potential The external potential properties.
 */
static INLINE void potential_print_backend(
    const struct external_potential* potential) {

  message("External potential is 'Turbulent driving'.");
  message("%i modes, dt = %g", potential->number_of_modes, potential->dt);
}

static INLINE void potential_update(double time,
                                    struct external_potential* potential) {

  /* first, check if we need to do anything */
  if (potential->number_of_steps * potential->dt < time) {

    /* reset the amplitudes */
    for (int i = 0; i < 6 * potential->number_of_modes; ++i) {
      potential->amplitudes[i] = 0.;
    }

    /* accumulate contributions to the forcing until we reach the desired point
       in time */
    while (potential->number_of_steps * potential->dt < time) {
      for (int i = 0; i < potential->number_of_modes; ++i) {

        const double phi =
            2. * M_PI * gsl_rng_uniform(potential->random_generator);
        const double theta1 =
            2. * M_PI * gsl_rng_uniform(potential->random_generator);
        const double theta2 =
            2. * M_PI * gsl_rng_uniform(potential->random_generator);

        const double ga = sin(phi);
        const double gb = cos(phi);
        const double real_rand1 = cos(theta1) * ga;
        const double imag_rand1 = sin(theta1) * ga;
        const double real_rand2 = cos(theta2) * gb;
        const double imag_rand2 = sin(theta2) * gb;

        const double kforce = potential->forcing[i];
        potential->amplitudes[6 * i + 0] =
            kforce * (real_rand1 * potential->unit_vectors[6 * i + 0] +
                      real_rand2 * potential->unit_vectors[6 * i + 3]);
        potential->amplitudes[6 * i + 1] =
            kforce * (real_rand1 * potential->unit_vectors[6 * i + 1] +
                      real_rand2 * potential->unit_vectors[6 * i + 4]);
        potential->amplitudes[6 * i + 2] =
            kforce * (real_rand1 * potential->unit_vectors[6 * i + 2] +
                      real_rand2 * potential->unit_vectors[6 * i + 5]);
        potential->amplitudes[6 * i + 3] =
            kforce * (imag_rand1 * potential->unit_vectors[6 * i + 0] +
                      imag_rand2 * potential->unit_vectors[6 * i + 3]);
        potential->amplitudes[6 * i + 4] =
            kforce * (imag_rand1 * potential->unit_vectors[6 * i + 1] +
                      imag_rand2 * potential->unit_vectors[6 * i + 4]);
        potential->amplitudes[6 * i + 5] =
            kforce * (imag_rand1 * potential->unit_vectors[6 * i + 2] +
                      imag_rand2 * potential->unit_vectors[6 * i + 5]);
      }
      ++potential->number_of_steps;
    }
  }
}

#endif /* SWIFT_POTENTIAL_TURBULENT_DRIVING_H */
