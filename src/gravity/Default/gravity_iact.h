/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#ifndef SWIFT_DEFAULT_GRAVITY_IACT_H
#define SWIFT_DEFAULT_GRAVITY_IACT_H

/* Includes. */
#include "const.h"
#include "kernel_gravity.h"
#include "kernel_long_gravity.h"
#include "multipole.h"
#include "vector.h"

/**
 * @brief Gravity forces between particles truncated by the long-range kernel
 */
__attribute__((always_inline)) INLINE static void runner_iact_grav_pp_truncated(
    float r2, const float *dx, struct gpart *gpi, struct gpart *gpj,
    float rlr_inv) {

  /* Apply the gravitational acceleration. */
  const float r = sqrtf(r2);
  const float ir = 1.f / r;
  const float mi = gpi->mass;
  const float mj = gpj->mass;
  const float hi = gpi->epsilon;
  const float hj = gpj->epsilon;
  const float u_lr = r * rlr_inv;
  float f_lr, fi, fj, W;

#ifdef SWIFT_DEBUG_CHECKS
  if (r == 0.f) error("Interacting particles with 0 distance");
#endif

  /* Get long-range correction */
  kernel_long_grav_eval(u_lr, &f_lr);

  if (r >= hi) {

    /* Get Newtonian gravity */
    fi = mj * ir * ir * ir * f_lr;

  } else {

    const float hi_inv = 1.f / hi;
    const float hi_inv3 = hi_inv * hi_inv * hi_inv;
    const float ui = r * hi_inv;

    kernel_grav_eval(ui, &W);

    /* Get softened gravity */
    fi = mj * hi_inv3 * W * f_lr;
  }

  if (r >= hj) {

    /* Get Newtonian gravity */
    fj = mi * ir * ir * ir * f_lr;

  } else {

    const float hj_inv = 1.f / hj;
    const float hj_inv3 = hj_inv * hj_inv * hj_inv;
    const float uj = r * hj_inv;

    kernel_grav_eval(uj, &W);

    /* Get softened gravity */
    fj = mi * hj_inv3 * W * f_lr;
  }

  const float fidx[3] = {fi * dx[0], fi * dx[1], fi * dx[2]};
  gpi->a_grav[0] -= fidx[0];
  gpi->a_grav[1] -= fidx[1];
  gpi->a_grav[2] -= fidx[2];

  const float fjdx[3] = {fj * dx[0], fj * dx[1], fj * dx[2]};
  gpj->a_grav[0] += fjdx[0];
  gpj->a_grav[1] += fjdx[1];
  gpj->a_grav[2] += fjdx[2];
}

/**
 * @brief Gravity forces between particles
 */
__attribute__((always_inline)) INLINE static void runner_iact_grav_pp(
    float r2, const float *dx, struct gpart *gpi, struct gpart *gpj) {

  /* Apply the gravitational acceleration. */
  const float r = sqrtf(r2);
  const float ir = 1.f / r;
  const float mi = gpi->mass;
  const float mj = gpj->mass;
  const float hi = gpi->epsilon;
  const float hj = gpj->epsilon;
  float fi, fj, W;

#ifdef SWIFT_DEBUG_CHECKS
  if (r == 0.f) error("Interacting particles with 0 distance");
#endif

  if (r >= hi) {

    /* Get Newtonian gravity */
    fi = mj * ir * ir * ir;

  } else {

    const float hi_inv = 1.f / hi;
    const float hi_inv3 = hi_inv * hi_inv * hi_inv;
    const float ui = r * hi_inv;

    kernel_grav_eval(ui, &W);

    /* Get softened gravity */
    fi = mj * hi_inv3 * W;
  }

  if (r >= hj) {

    /* Get Newtonian gravity */
    fj = mi * ir * ir * ir;

  } else {

    const float hj_inv = 1.f / hj;
    const float hj_inv3 = hj_inv * hj_inv * hj_inv;
    const float uj = r * hj_inv;

    kernel_grav_eval(uj, &W);

    /* Get softened gravity */
    fj = mi * hj_inv3 * W;
  }

  const float fidx[3] = {fi * dx[0], fi * dx[1], fi * dx[2]};
  gpi->a_grav[0] -= fidx[0];
  gpi->a_grav[1] -= fidx[1];
  gpi->a_grav[2] -= fidx[2];

  const float fjdx[3] = {fj * dx[0], fj * dx[1], fj * dx[2]};
  gpj->a_grav[0] += fjdx[0];
  gpj->a_grav[1] += fjdx[1];
  gpj->a_grav[2] += fjdx[2];
}

/**
 * @brief Gravity forces between particles truncated by the long-range kernel
 * (non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void
runner_iact_grav_pp_truncated_nonsym(float r2, const float *dx,
                                     struct gpart *gpi, const struct gpart *gpj,
                                     float rlr_inv) {

  /* Apply the gravitational acceleration. */
  const float r = sqrtf(r2);
  const float ir = 1.f / r;
  const float mj = gpj->mass;
  const float hi = gpi->epsilon;
  const float u_lr = r * rlr_inv;
  float f_lr, f, W;

#ifdef SWIFT_DEBUG_CHECKS
  if (r == 0.f) error("Interacting particles with 0 distance");
#endif

  /* Get long-range correction */
  kernel_long_grav_eval(u_lr, &f_lr);

  if (r >= hi) {

    /* Get Newtonian gravity */
    f = mj * ir * ir * ir * f_lr;

  } else {

    const float hi_inv = 1.f / hi;
    const float hi_inv3 = hi_inv * hi_inv * hi_inv;
    const float ui = r * hi_inv;

    kernel_grav_eval(ui, &W);

    /* Get softened gravity */
    f = mj * hi_inv3 * W * f_lr;
  }

  const float fdx[3] = {f * dx[0], f * dx[1], f * dx[2]};

  gpi->a_grav[0] -= fdx[0];
  gpi->a_grav[1] -= fdx[1];
  gpi->a_grav[2] -= fdx[2];
}

/**
 * @brief Gravity forces between particles (non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_grav_pp_nonsym(
    float r2, const float *dx, struct gpart *gpi, const struct gpart *gpj) {

  /* Apply the gravitational acceleration. */
  const float r = sqrtf(r2);
  const float ir = 1.f / r;
  const float mj = gpj->mass;
  const float hi = gpi->epsilon;
  float f, W;

#ifdef SWIFT_DEBUG_CHECKS
  if (r == 0.f) error("Interacting particles with 0 distance");
#endif

  if (r >= hi) {

    /* Get Newtonian gravity */
    f = mj * ir * ir * ir;

  } else {

    const float hi_inv = 1.f / hi;
    const float hi_inv3 = hi_inv * hi_inv * hi_inv;
    const float ui = r * hi_inv;

    kernel_grav_eval(ui, &W);

    /* Get softened gravity */
    f = mj * hi_inv3 * W;
  }

  const float fdx[3] = {f * dx[0], f * dx[1], f * dx[2]};

  gpi->a_grav[0] -= fdx[0];
  gpi->a_grav[1] -= fdx[1];
  gpi->a_grav[2] -= fdx[2];
}

#endif /* SWIFT_DEFAULT_GRAVITY_IACT_H */
