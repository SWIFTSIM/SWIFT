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

/* Config parameters. */
#include "turbulence.h"

#include "../config.h"
#include "engine.h"

#ifdef TURBULENCE_DRIVING_ALVELIUS

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
void turbulence_init_backend(struct swift_params* parameter_file,
                             const struct phys_const* phys_const,
                             const struct unit_system* us,
                             const struct space* s,
                             struct turbulence_driving* turbulence) {

  /* make sure the box is a cube */
  if (s->dim[0] != s->dim[1] || s->dim[0] != s->dim[2]) {
    error("Turbulent forcing only works in a cubic box!");
  }

  /* get dimensionless parameters */
  const int random_seed = parser_get_opt_param_int(
      parameter_file, "TurbulentDriving:random_seed", 42);
  const double kmin =
      parser_get_opt_param_double(parameter_file, "TurbulentDriving:kmin", 2.);
  const double kmax =
      parser_get_opt_param_double(parameter_file, "TurbulentDriving:kmax", 3.);
  const double kforcing = parser_get_opt_param_double(
      parameter_file, "TurbulentDriving:kforcing", 2.5);
  const double concentration_factor = parser_get_opt_param_double(
      parameter_file, "TurbulentDriving:concentration_factor", 0.2);

  /* get parameters with units */
  const double power_forcing_cgs = parser_get_opt_param_double(
      parameter_file, "TurbulentDriving:power_forcing_in_cm2_per_s3", 17.);
  const double dtfor_cgs = parser_get_opt_param_double(
      parameter_file, "TurbulentDriving:dt_forcing_in_s", 1.e6);
  const double starting_time_cgs = parser_get_opt_param_double(
      parameter_file, "TurbulentDriving:starting_time_in_s", 0.);

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
  turbulence->random_generator = gsl_rng_alloc(gsl_rng_ranlxd2);
  gsl_rng_set(turbulence->random_generator, random_seed);

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
  turbulence->k =
      (double*)swift_malloc("turbulent_k", 3 * num_modes * sizeof(double));
  turbulence->amplitudes =
      (double*)swift_malloc("turbulent_A", 6 * num_modes * sizeof(double));
  turbulence->unit_vectors =
      (double*)swift_malloc("turbulent_uv", 6 * num_modes * sizeof(double));
  turbulence->forcing =
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
          const double kdiff = k - kforcing;
          const double sqrtk12 = sqrt(pwrk1 + pwrk2);
          const double invkk = 1. / kk;
          const double invk = 1. / k;
          if (sqrtk12 > 0.) {
            const double invsqrtk12 = 1. / sqrtk12;
            turbulence->unit_vectors[6 * kindex + 0] = k2 * invsqrtk12;
            turbulence->unit_vectors[6 * kindex + 1] = -k1 * invsqrtk12;
            turbulence->unit_vectors[6 * kindex + 2] = 0.;
            turbulence->unit_vectors[6 * kindex + 3] =
                k1 * k3 * invsqrtk12 * invk;
            turbulence->unit_vectors[6 * kindex + 4] =
                k2 * k3 * invsqrtk12 * invk;
            turbulence->unit_vectors[6 * kindex + 5] = -sqrtk12 * invk;
          } else {
            const double sqrtk13 = sqrt(pwrk1 + pwrk3);
            const double invsqrtk13 = 1. / sqrtk13;
            turbulence->unit_vectors[6 * kindex + 0] = -k3 * invsqrtk13;
            turbulence->unit_vectors[6 * kindex + 1] = 0.;
            turbulence->unit_vectors[6 * kindex + 2] = k1 * invsqrtk13;
            turbulence->unit_vectors[6 * kindex + 3] =
                k1 * k2 * invsqrtk13 * invk;
            turbulence->unit_vectors[6 * kindex + 4] = -sqrtk13 * invk;
            turbulence->unit_vectors[6 * kindex + 5] =
                k2 * k3 * invsqrtk13 * invk;
          }

          turbulence->k[3 * kindex + 0] = k1 * Linv;
          turbulence->k[3 * kindex + 1] = k2 * Linv;
          turbulence->k[3 * kindex + 2] = k3 * Linv;
          const double gaussian_spectrum = exp(-kdiff * kdiff * cinv) * invkk;
          spectrum_sum += gaussian_spectrum;
          turbulence->forcing[kindex] = gaussian_spectrum;
          ++kindex;
        }
      }
    }
  }

  /* normalise the forcing */
  const double norm = power_forcing / (spectrum_sum * dtfor);
  for (int i = 0; i < num_modes; ++i) {
    turbulence->forcing[i] *= norm;
    turbulence->forcing[i] = sqrt(turbulence->forcing[i]);
  }

  /* fast-forward the driving to the desired point in time */
  int num_steps = 0;
  while (num_steps * dtfor < starting_time) {
    /* 3 random numbers are generated per mode in potential_update() */
    for (int i = 0; i < 3 * num_modes; ++i) {
      gsl_rng_uniform(turbulence->random_generator);
    }
    ++num_steps;
  }

  /* initialise the final variables */
  turbulence->number_of_steps = 0;
  turbulence->dt = dtfor;
  turbulence->number_of_modes = num_modes;
}

/**
 * @brief Prints the properties of the external potential to stdout.
 *
 * @param  potential The external potential properties.
 */
void turbulence_print_backend(const struct turbulence_driving* turbulence) {

  message("Turbulence driving mode is 'Alvelius'.");
  message("%i modes, dt = %g", turbulence->number_of_modes, turbulence->dt);
}

void turbulence_update(struct engine* restrict e) {

  const double time = e->time;
  struct turbulence_driving* turbulence = e->turbulence;

  /* first, check if we need to do anything */
  if (turbulence->number_of_steps * turbulence->dt < time) {

    /* reset the amplitudes */
    for (int i = 0; i < 6 * turbulence->number_of_modes; ++i) {
      turbulence->amplitudes[i] = 0.;
    }

    /* accumulate contributions to the forcing until we reach the desired point
       in time */
    while (turbulence->number_of_steps * turbulence->dt < time) {
      for (int i = 0; i < turbulence->number_of_modes; ++i) {

        const double phi =
            2. * M_PI * gsl_rng_uniform(turbulence->random_generator);
        const double theta1 =
            2. * M_PI * gsl_rng_uniform(turbulence->random_generator);
        const double theta2 =
            2. * M_PI * gsl_rng_uniform(turbulence->random_generator);

        const double ga = sin(phi);
        const double gb = cos(phi);
        const double real_rand1 = cos(theta1) * ga;
        const double imag_rand1 = sin(theta1) * ga;
        const double real_rand2 = cos(theta2) * gb;
        const double imag_rand2 = sin(theta2) * gb;

        const double kforce = turbulence->forcing[i];
        turbulence->amplitudes[6 * i + 0] +=
            kforce * (real_rand1 * turbulence->unit_vectors[6 * i + 0] +
                      real_rand2 * turbulence->unit_vectors[6 * i + 3]);
        turbulence->amplitudes[6 * i + 1] +=
            kforce * (real_rand1 * turbulence->unit_vectors[6 * i + 1] +
                      real_rand2 * turbulence->unit_vectors[6 * i + 4]);
        turbulence->amplitudes[6 * i + 2] +=
            kforce * (real_rand1 * turbulence->unit_vectors[6 * i + 2] +
                      real_rand2 * turbulence->unit_vectors[6 * i + 5]);
        turbulence->amplitudes[6 * i + 3] +=
            kforce * (imag_rand1 * turbulence->unit_vectors[6 * i + 0] +
                      imag_rand2 * turbulence->unit_vectors[6 * i + 3]);
        turbulence->amplitudes[6 * i + 4] +=
            kforce * (imag_rand1 * turbulence->unit_vectors[6 * i + 1] +
                      imag_rand2 * turbulence->unit_vectors[6 * i + 4]);
        turbulence->amplitudes[6 * i + 5] +=
            kforce * (imag_rand1 * turbulence->unit_vectors[6 * i + 2] +
                      imag_rand2 * turbulence->unit_vectors[6 * i + 5]);
      }
      ++turbulence->number_of_steps;
    }
  }

  const int count = e->s->nr_parts;
  struct part* parts = e->s->parts;
  struct xpart* xparts = e->s->xparts;
  for (int i = 0; i < count; ++i) {
    struct part* p = &parts[i];
    struct xpart* xp = &xparts[i];
    turbulence_accelerate(p, xp, turbulence);
  }
}

void turbulence_accelerate(struct part* restrict p, struct xpart* restrict xp,
                           const struct turbulence_driving* restrict
                               turbulence) {

  const double* x = p->x;

  double force[3] = {0., 0., 0.};
  for (int ik = 0; ik < turbulence->number_of_modes; ++ik) {
    const double fr[3] = {turbulence->amplitudes[6 * ik + 0],
                          turbulence->amplitudes[6 * ik + 1],
                          turbulence->amplitudes[6 * ik + 2]};
    const double fi[3] = {turbulence->amplitudes[6 * ik + 3],
                          turbulence->amplitudes[6 * ik + 4],
                          turbulence->amplitudes[6 * ik + 5]};
    const double k[3] = {turbulence->k[3 * ik + 0], turbulence->k[3 * ik + 1],
                         turbulence->k[3 * ik + 2]};

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

  p->v[0] += force[0] * turbulence->dt;
  p->v[1] += force[1] * turbulence->dt;
  p->v[2] += force[2] * turbulence->dt;
  xp->v_full[0] += force[0] * turbulence->dt;
  xp->v_full[1] += force[1] * turbulence->dt;
  xp->v_full[2] += force[2] * turbulence->dt;
}

#else /* TURBULENCE_DRIVING_NONE */

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
void turbulence_init_backend(struct swift_params* parameter_file,
                             const struct phys_const* phys_const,
                             const struct unit_system* us,
                             const struct space* s,
                             struct turbulence_driving* turbulence) {}

/**
 * @brief Prints the properties of the external potential to stdout.
 *
 * @param  potential The external potential properties.
 */
void turbulence_print_backend(const struct turbulence_driving* turbulence) {}

void turbulence_update(struct engine* restrict e) {}

void turbulence_accelerate(struct part* restrict p, struct xpart* restrict xp,
                           const struct turbulence_driving* restrict
                               turbulence) {}

#endif /* TURBULENCE_DRIVING */
