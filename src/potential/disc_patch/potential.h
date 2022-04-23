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
#ifndef SWIFT_DISC_PATCH_H
#define SWIFT_DISC_PATCH_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <math.h>

/* Local includes. */
#include "error.h"
#include "gravity.h"
#include "minmax.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"

/**
 * @brief External Potential Properties - Disc patch case
 *
 * See Creasey, Theuns & Bower, 2013, MNRAS, Volume 429, Issue 3, p.1922-1948.
 *
 * We truncate the accelerations beyond z_trunc using a 1-cos(z) function
 * that smoothly brings the accelerations to 0 at z_max.
 */
struct external_potential {

  /*! Surface density of the disc (sigma) */
  float surface_density;

  /*! Disc scale-height (b) */
  float scale_height;

  /*! Inverse of disc scale-height (1/b) */
  float scale_height_inv;

  /*! Position of the disc along the x-axis */
  float x_disc;

  /*! Position above which the accelerations get truncated */
  float x_trunc;

  /*! Position above which the accelerations are zero */
  float x_max;

  /*! The truncated transition regime */
  float x_trans;

  /*! Inverse of the truncated transition regime */
  float x_trans_inv;

  /*! Dynamical time of the system */
  float dynamical_time;

  /*! Time over which to grow the disk */
  float growth_time;

  /*! Inverse of the growth time */
  float growth_time_inv;

  /*! Time-step condition pre-factor */
  float timestep_mult;

  /*! Constant pre-factor (2 pi G sigma) */
  float norm;

  /*! Constant pre-factor (2 pi sigma)*/
  float norm_over_G;
};

/**
 * @brief Computes the time-step from the acceleration due to a hydrostatic
 * disc.
 *
 * See Creasey, Theuns & Bower, 2013, MNRAS, Volume 429, Issue 3, p.1922-1948,
 * equations 17 and 20.
 * We do not use the truncated potential here.
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

  /* initilize time step to disc dynamical time */
  const float dt_dyn = potential->dynamical_time;
  const float b = potential->scale_height;
  const float b_inv = potential->scale_height_inv;
  const float norm = potential->norm;

  /* absolute value of height above disc */
  const float dx = fabs(g->x[0] - potential->x_disc);

  /* vertical acceleration */
  const float x_accel = norm * tanhf(dx * b_inv);

  float dt = dt_dyn;

  /* demand that dt * velocity <  fraction of scale height of disc */
  if (g->v_full[0] != 0.f) {

    const float dt1 = b / fabsf(g->v_full[0]);
    dt = min(dt1, dt);
  }

  /* demand that dt^2 * acceleration < fraction of scale height of disc */
  if (x_accel != 0.f) {

    const float dt2 = b / fabsf(x_accel);
    if (dt2 < dt * dt) dt = sqrtf(dt2);
  }

  /* demand that dt^3 * jerk < fraction of scale height of disc */
  if (g->v_full[0] != 0.f) {

    const float cosh_dx_inv = 1.f / coshf(dx * b_inv);
    const float cosh_dx_inv2 = cosh_dx_inv * cosh_dx_inv;
    const float dx_accel_over_dt =
        norm * cosh_dx_inv2 * b_inv * fabsf(g->v_full[0]);

    const float dt3 = b / fabsf(dx_accel_over_dt);
    if (dt3 < dt * dt * dt) dt = cbrtf(dt3);
  }

  return potential->timestep_mult * dt;
}

/**
 * @brief Computes the gravitational acceleration along x due to a hydrostatic
 * disc
 *
 * See Creasey, Theuns & Bower, 2013, MNRAS, Volume 429, Issue 3, p.1922-1948,
 * equation 17.
 * We truncate the accelerations beyond x_trunc using a 1-cos(x) function
 * that smoothly brings the accelerations to 0 at x_max.
 *
 * @param time The current time in internal units.
 * @param potential The properties of the potential.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static void external_gravity_acceleration(
    double time, const struct external_potential* restrict potential,
    const struct phys_const* restrict phys_const, struct gpart* restrict g) {

  const float dx = g->x[0] - potential->x_disc;
  const float abs_dx = fabsf(dx);
  const float t_growth = potential->growth_time;
  const float t_growth_inv = potential->growth_time_inv;
  const float b = potential->scale_height;
  const float b_inv = potential->scale_height_inv;
  const float x_trunc = potential->x_trunc;
  const float x_max = potential->x_max;
  const float x_trans_inv = potential->x_trans_inv;
  const float norm_over_G = potential->norm_over_G;

  /* Are we still growing the disc ? */
  const float reduction_factor = time < t_growth ? time * t_growth_inv : 1.f;

  /* Truncated or not ? */
  float a_x;
  float pot;
  if (abs_dx < x_trunc) {

    /* Acc. 2 pi sigma tanh(x/b) */
    a_x = reduction_factor * norm_over_G * tanhf(abs_dx * b_inv);
    pot = -reduction_factor * norm_over_G * logf(coshf(abs_dx * b_inv)) * b;
  } else if (abs_dx < x_max) {

    /* Acc. 2 pi sigma tanh(x/b) [1/2 + 1/2cos((x-xmax)/(pi x_trans))] */
    a_x =
        reduction_factor * norm_over_G * tanhf(abs_dx * b_inv) *
        (0.5f + 0.5f * cosf((float)(M_PI) * (abs_dx - x_trunc) * x_trans_inv));
    pot = 0.f;
  } else {

    /* Acc. 0 */
    a_x = 0.f;
    pot = 0.f;
  }

  /* Get the correct sign. Recall G is multipiled in later on */
  if (dx > 0) g->a_grav[0] -= a_x;
  if (dx < 0) g->a_grav[0] += a_x;
  gravity_add_comoving_potential(g, pot);
}

/**
 * @brief Computes the gravitational potential energy of a particle in the
 * disc patch potential.
 *
 * See Creasey, Theuns & Bower, 2013, MNRAS, Volume 429, Issue 3, p.1922-1948,
 * equation 22.
 * We truncate the accelerations beyond x_trunc using a 1-cos(x) function
 * that smoothly brings the accelerations to 0 at x_max.
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

  const float dx = gp->x[0] - potential->x_disc;
  const float abs_dx = fabsf(dx);
  const float t_growth = potential->growth_time;
  const float t_growth_inv = potential->growth_time_inv;
  const float b = potential->scale_height;
  const float b_inv = potential->scale_height_inv;
  const float norm = potential->norm;
  const float x_trunc = potential->x_trunc;
  const float x_max = potential->x_max;

  /* Are we still growing the disc ? */
  const float reduction_factor = time < t_growth ? time * t_growth_inv : 1.f;

  /* Truncated or not ? */
  float pot;
  if (abs_dx < x_trunc) {

    /* Potential (2 pi G sigma b ln(cosh(x/b)) */
    pot = b * logf(coshf(dx * b_inv));
  } else if (abs_dx < x_max) {

    /* Potential. At x>>b, phi(x) = norm * x / b */
    pot = 0.f;

  } else {

    pot = 0.f;
  }

  return pot * reduction_factor * norm;
}

/**
 * @brief Initialises the external potential properties in the internal system
 * of units.
 *
 * See Creasey, Theuns & Bower, 2013, MNRAS, Volume 429, Issue 3, p.1922-1948,
 * equation 22.
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

  potential->surface_density = parser_get_param_double(
      parameter_file, "DiscPatchPotential:surface_density");
  potential->scale_height = parser_get_param_double(
      parameter_file, "DiscPatchPotential:scale_height");
  potential->x_disc =
      parser_get_param_double(parameter_file, "DiscPatchPotential:x_disc");
  potential->x_trunc = parser_get_opt_param_double(
      parameter_file, "DiscPatchPotential:x_trunc", FLT_MAX);
  potential->x_max = parser_get_opt_param_double(
      parameter_file, "DiscPatchPotential:x_max", FLT_MAX);
  potential->timestep_mult = parser_get_param_double(
      parameter_file, "DiscPatchPotential:timestep_mult");
  potential->growth_time = parser_get_opt_param_double(
      parameter_file, "DiscPatchPotential:growth_time", 0.);

  /* Compute the dynamical time */
  potential->dynamical_time =
      sqrt(potential->scale_height /
           (phys_const->const_newton_G * potential->surface_density));

  /* Convert the growth time multiplier to physical time */
  potential->growth_time *= potential->dynamical_time;

  /* Some cross-checks */
  if (potential->x_trunc > potential->x_max)
    error("Potential truncation x larger than maximal z");
  if (potential->x_trunc < potential->scale_height)
    error("Potential truncation x smaller than scale height");

  /* Compute derived quantities */
  potential->scale_height_inv = 1. / potential->scale_height;
  potential->norm =
      2. * M_PI * phys_const->const_newton_G * potential->surface_density;
  potential->norm_over_G = 2 * M_PI * potential->surface_density;
  potential->x_trans = potential->x_max - potential->x_trunc;

  if (potential->x_trans != 0.f)
    potential->x_trans_inv = 1. / potential->x_trans;
  else
    potential->x_trans_inv = FLT_MAX;

  if (potential->growth_time != 0.)
    potential->growth_time_inv = 1. / potential->growth_time;
  else
    potential->growth_time_inv = FLT_MAX;
}

/**
 * @brief Prints the properties of the external potential to stdout.
 *
 * @param  potential The external potential properties.
 */
static INLINE void potential_print_backend(
    const struct external_potential* potential) {

  message(
      "External potential is 'Disk-patch' with Sigma=%f, x_disc=%f, b=%f and "
      "dt_mult=%f.",
      potential->surface_density, potential->x_disc, potential->scale_height,
      potential->timestep_mult);

  if (potential->x_max < FLT_MAX)
    message("Potential will be truncated at x_trunc=%f and zeroed at x_max=%f",
            potential->x_trunc, potential->x_max);

  if (potential->growth_time > 0.)
    message("Disc will grow for %f [time_units]. (%f dynamical time)",
            potential->growth_time,
            potential->growth_time / potential->dynamical_time);
}

#endif /* SWIFT_DISC_PATCH_H */
