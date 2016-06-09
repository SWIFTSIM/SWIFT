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
#ifndef SWIFT_RUNNER_IACT_GRAV_H
#define SWIFT_RUNNER_IACT_GRAV_H

/* Includes. */
#include "const.h"
#include "kernel_gravity.h"
#include "multipole.h"
#include "vector.h"

/**
 * @brief Gravity forces between particles
 */
__attribute__((always_inline)) INLINE static void runner_iact_grav_pp(
    float r2, const float *dx, struct gpart *gpi, struct gpart *gpj) {

  /* Apply the gravitational acceleration. */
  const float ir = 1.0f / sqrtf(r2);
  const float w = ir * ir * ir;
  const float wdx[3] = {w * dx[0], w * dx[1], w * dx[2]};
  const float mi = gpi->mass;
  const float mj = gpj->mass;

  gpi->a_grav[0] -= wdx[0] * mj;
  gpi->a_grav[1] -= wdx[1] * mj;
  gpi->a_grav[2] -= wdx[2] * mj;
  gpi->mass_interacted += mj;

  gpj->a_grav[0] += wdx[0] * mi;
  gpj->a_grav[1] += wdx[1] * mi;
  gpj->a_grav[2] += wdx[2] * mi;
  gpj->mass_interacted += mi;
}

__attribute__((always_inline)) INLINE static void runner_iact_grav_pp_nonsym(
    float r2, const float *dx, struct gpart *gpi, const struct gpart *gpj) {

  /* Apply the gravitational acceleration. */
  const float ir = 1.0f / sqrtf(r2);
  const float w = ir * ir * ir;
  const float wdx[3] = {w * dx[0], w * dx[1], w * dx[2]};
  const float mj = gpj->mass;

  gpi->a_grav[0] -= wdx[0] * mj;
  gpi->a_grav[1] -= wdx[1] * mj;
  gpi->a_grav[2] -= wdx[2] * mj;
  gpi->mass_interacted += mj;
}

/**
 * @brief Gravity forces between particles
 */
__attribute__((always_inline)) INLINE static void runner_iact_grav_pm(
    float r2, const float *dx, struct gpart *gp,
    const struct multipole *multi) {

  /* Apply the gravitational acceleration. */
  const float r = sqrtf(r2);
  const float ir = 1.f / r;
  const float mrinv3 = multi->mass * ir * ir * ir;

#if multipole_order < 2

  /* 0th and 1st order terms */
  gp->a_grav[0] += mrinv3 * dx[0];
  gp->a_grav[1] += mrinv3 * dx[1];
  gp->a_grav[2] += mrinv3 * dx[2];

  gp->mass_interacted += multi->mass;
#elif multipole_order == 2
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

  gp->a_grav[0] -= C * dx[0] + D2 * qRx;
  gp->a_grav[1] -= C * dx[1] + D2 * qRy;
  gp->a_grav[2] -= C * dx[2] + D2 * qRz;

  gp->mass_interacted += multi->mass;
#else
#error "Multipoles of order >2 not yet implemented."
#endif
}

#endif /* SWIFT_RUNNER_IACT_GRAV_H */
