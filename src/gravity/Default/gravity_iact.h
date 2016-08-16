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
 * @brief Gravity forces between particles
 */
__attribute__((always_inline)) INLINE static void runner_iact_grav_pp(
    float rlr_inv, float r2, const float *dx, struct gpart *gpi,
    struct gpart *gpj) {

  /* Apply the gravitational acceleration. */
  const float r = sqrtf(r2);
  const float ir = 1.f / r;
  const float mi = gpi->mass;
  const float mj = gpj->mass;
  const float hi = gpi->epsilon;
  const float hj = gpj->epsilon;
  const float u = r * rlr_inv;
  float f_lr, fi, fj, W;

  /* Get long-range correction */
  kernel_long_grav_eval(u, &f_lr);

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
 * @brief Gravity forces between particles (non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_grav_pp_nonsym(
    float rlr_inv, float r2, const float *dx, struct gpart *gpi,
    const struct gpart *gpj) {

  /* Apply the gravitational acceleration. */
  const float r = sqrtf(r2);
  const float ir = 1.f / r;
  const float mj = gpj->mass;
  const float hi = gpi->epsilon;
  const float u = r * rlr_inv;
  float f_lr, f, W;

  /* Get long-range correction */
  kernel_long_grav_eval(u, &f_lr);

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
 * @brief Gravity forces between particle and multipole
 */
__attribute__((always_inline)) INLINE static void runner_iact_grav_pm(
    float rlr_inv, float r2, const float *dx, struct gpart *gp,
    const struct multipole *multi) {

  /* Apply the gravitational acceleration. */
  const float r = sqrtf(r2);
  const float ir = 1.f / r;
  const float u = r * rlr_inv;
  const float mrinv3 = multi->mass * ir * ir * ir;

  /* Get long-range correction */
  float f_lr;
  kernel_long_grav_eval(u, &f_lr);

#if const_gravity_multipole_order < 2

  /* 0th and 1st order terms */
  gp->a_grav[0] += mrinv3 * f_lr * dx[0];
  gp->a_grav[1] += mrinv3 * f_lr * dx[1];
  gp->a_grav[2] += mrinv3 * f_lr * dx[2];

#elif const_gravity_multipole_order == 2
  /* Terms up to 2nd order (quadrupole) */

  /* Follows the notation in Bonsai */
  const float mrinv5 = mrinv3 * ir * ir;
  const float mrinv7 = mrinv5 * ir * ir;

  const float D1 = -mrinv3;
  const float D2 = 3.f * mrinv5;
  const float D3 = -15.f * mrinv7;

  const float q = multi->I_xx + multi->I_yy + multi->I_zz;
  const float qRx =
      multi->I_xx * dx[0] + multi->I_xy * dx[1] + multi->I_xz * dx[2];
  const float qRy =
      multi->I_xy * dx[0] + multi->I_yy * dx[1] + multi->I_yz * dx[2];
  const float qRz =
      multi->I_xz * dx[0] + multi->I_yz * dx[1] + multi->I_zz * dx[2];
  const float qRR = qRx * dx[0] + qRy * dx[1] + qRz * dx[2];
  const float C = D1 + 0.5f * D2 * q + 0.5f * D3 * qRR;

  gp->a_grav[0] -= f_lr * (C * dx[0] + D2 * qRx);
  gp->a_grav[1] -= f_lr * (C * dx[1] + D2 * qRy);
  gp->a_grav[2] -= f_lr * (C * dx[2] + D2 * qRz);

#else
#error "Multipoles of order >2 not yet implemented."
#endif
}

#endif /* SWIFT_DEFAULT_GRAVITY_IACT_H */
