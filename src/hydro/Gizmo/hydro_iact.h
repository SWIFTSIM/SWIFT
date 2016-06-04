/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *                    Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
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

#include "riemann.h"

#define USE_GRADIENTS
#define PER_FACE_LIMITER
/* #define PRINT_ID 0 */

/* this corresponds to task_subtype_hydro_loop1 */
__attribute__((always_inline)) INLINE static void runner_iact_hydro_loop1(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float r = sqrtf(r2);
  float xi, xj;
  float h_inv;
  float wi, wj, wi_dx, wj_dx;
  int k, l;

  /* Compute density of pi. */
  h_inv = 1.0 / hi;
  xi = r * h_inv;
  kernel_deval(xi, &wi, &wi_dx);

  pi->density.wcount += wi;
  pi->density.wcount_dh -= xi * wi_dx;

  /* these are eqns. (1) and (2) in the summary */
  pi->geometry.volume += wi;
  for (k = 0; k < 3; k++)
    for (l = 0; l < 3; l++) pi->geometry.matrix_E[k][l] += dx[k] * dx[l] * wi;

#ifdef SPH_GRADIENTS
  /* very basic gradient estimate */
  pi->primitives.gradients.rho[0] -=
      wi_dx * dx[0] * (pi->primitives.rho - pj->primitives.rho) / r;
  pi->primitives.gradients.rho[1] -=
      wi_dx * dx[1] * (pi->primitives.rho - pj->primitives.rho) / r;
  pi->primitives.gradients.rho[2] -=
      wi_dx * dx[2] * (pi->primitives.rho - pj->primitives.rho) / r;

  pi->primitives.gradients.v[0][0] -=
      wi_dx * dx[0] * (pi->primitives.v[0] - pj->primitives.v[0]) / r;
  pi->primitives.gradients.v[0][1] -=
      wi_dx * dx[1] * (pi->primitives.v[0] - pj->primitives.v[0]) / r;
  pi->primitives.gradients.v[0][2] -=
      wi_dx * dx[2] * (pi->primitives.v[0] - pj->primitives.v[0]) / r;

  pi->primitives.gradients.v[1][0] -=
      wi_dx * dx[0] * (pi->primitives.v[1] - pj->primitives.v[1]) / r;
  pi->primitives.gradients.v[1][1] -=
      wi_dx * dx[1] * (pi->primitives.v[1] - pj->primitives.v[1]) / r;
  pi->primitives.gradients.v[1][2] -=
      wi_dx * dx[2] * (pi->primitives.v[1] - pj->primitives.v[1]) / r;

  pi->primitives.gradients.v[2][0] -=
      wi_dx * dx[0] * (pi->primitives.v[2] - pj->primitives.v[2]) / r;
  pi->primitives.gradients.v[2][1] -=
      wi_dx * dx[1] * (pi->primitives.v[2] - pj->primitives.v[2]) / r;
  pi->primitives.gradients.v[2][2] -=
      wi_dx * dx[2] * (pi->primitives.v[2] - pj->primitives.v[2]) / r;

  pi->primitives.gradients.P[0] -=
      wi_dx * dx[0] * (pi->primitives.P - pj->primitives.P) / r;
  pi->primitives.gradients.P[1] -=
      wi_dx * dx[1] * (pi->primitives.P - pj->primitives.P) / r;
  pi->primitives.gradients.P[2] -=
      wi_dx * dx[2] * (pi->primitives.P - pj->primitives.P) / r;

  /* basic slope limiter: collect the maximal and the minimal value for the
   * primitive variables among the ngbs */
  pi->primitives.limiter.rho[0] =
      fmin(pj->primitives.rho, pi->primitives.limiter.rho[0]);
  pi->primitives.limiter.rho[1] =
      fmax(pj->primitives.rho, pi->primitives.limiter.rho[1]);

  pi->primitives.limiter.v[0][0] =
      fmin(pj->primitives.v[0], pi->primitives.limiter.v[0][0]);
  pi->primitives.limiter.v[0][1] =
      fmax(pj->primitives.v[0], pi->primitives.limiter.v[0][1]);
  pi->primitives.limiter.v[1][0] =
      fmin(pj->primitives.v[1], pi->primitives.limiter.v[1][0]);
  pi->primitives.limiter.v[1][1] =
      fmax(pj->primitives.v[1], pi->primitives.limiter.v[1][1]);
  pi->primitives.limiter.v[2][0] =
      fmin(pj->primitives.v[2], pi->primitives.limiter.v[2][0]);
  pi->primitives.limiter.v[2][1] =
      fmax(pj->primitives.v[2], pi->primitives.limiter.v[2][1]);

  pi->primitives.limiter.P[0] =
      fmin(pj->primitives.P, pi->primitives.limiter.P[0]);
  pi->primitives.limiter.P[1] =
      fmax(pj->primitives.P, pi->primitives.limiter.P[1]);

  pi->primitives.limiter.maxr = fmax(r, pi->primitives.limiter.maxr);
#endif

  /* Compute density of pj. */
  h_inv = 1.0 / hj;
  xj = r * h_inv;
  kernel_deval(xj, &wj, &wj_dx);

  pj->density.wcount += wj;
  pj->density.wcount_dh -= xj * wj_dx;

  /* these are eqns. (1) and (2) in the summary */
  pj->geometry.volume += wj;
  for (k = 0; k < 3; k++)
    for (l = 0; l < 3; l++) pj->geometry.matrix_E[k][l] += dx[k] * dx[l] * wj;

#ifdef SPH_GRADIENTS
  /* very basic gradient estimate */
  /* signs are the same as before, since we swap i and j twice */
  pj->primitives.gradients.rho[0] -=
      wj_dx * dx[0] * (pi->primitives.rho - pj->primitives.rho) / r;
  pj->primitives.gradients.rho[1] -=
      wj_dx * dx[1] * (pi->primitives.rho - pj->primitives.rho) / r;
  pj->primitives.gradients.rho[2] -=
      wj_dx * dx[2] * (pi->primitives.rho - pj->primitives.rho) / r;

  pj->primitives.gradients.v[0][0] -=
      wj_dx * dx[0] * (pi->primitives.v[0] - pj->primitives.v[0]) / r;
  pj->primitives.gradients.v[0][1] -=
      wj_dx * dx[1] * (pi->primitives.v[0] - pj->primitives.v[0]) / r;
  pj->primitives.gradients.v[0][2] -=
      wj_dx * dx[2] * (pi->primitives.v[0] - pj->primitives.v[0]) / r;

  pj->primitives.gradients.v[1][0] -=
      wj_dx * dx[0] * (pi->primitives.v[1] - pj->primitives.v[1]) / r;
  pj->primitives.gradients.v[1][1] -=
      wj_dx * dx[1] * (pi->primitives.v[1] - pj->primitives.v[1]) / r;
  pj->primitives.gradients.v[1][2] -=
      wj_dx * dx[2] * (pi->primitives.v[1] - pj->primitives.v[1]) / r;

  pj->primitives.gradients.v[2][0] -=
      wj_dx * dx[0] * (pi->primitives.v[2] - pj->primitives.v[2]) / r;
  pj->primitives.gradients.v[2][1] -=
      wj_dx * dx[1] * (pi->primitives.v[2] - pj->primitives.v[2]) / r;
  pj->primitives.gradients.v[2][2] -=
      wj_dx * dx[2] * (pi->primitives.v[2] - pj->primitives.v[2]) / r;

  pj->primitives.gradients.P[0] -=
      wj_dx * dx[0] * (pi->primitives.P - pj->primitives.P) / r;
  pj->primitives.gradients.P[1] -=
      wj_dx * dx[1] * (pi->primitives.P - pj->primitives.P) / r;
  pj->primitives.gradients.P[2] -=
      wj_dx * dx[2] * (pi->primitives.P - pj->primitives.P) / r;

  /* basic slope limiter: collect the maximal and the minimal value for the
   * primitive variables among the ngbs */
  pj->primitives.limiter.rho[0] =
      fmin(pi->primitives.rho, pj->primitives.limiter.rho[0]);
  pj->primitives.limiter.rho[1] =
      fmax(pi->primitives.rho, pj->primitives.limiter.rho[1]);

  pj->primitives.limiter.v[0][0] =
      fmin(pi->primitives.v[0], pj->primitives.limiter.v[0][0]);
  pj->primitives.limiter.v[0][1] =
      fmax(pi->primitives.v[0], pj->primitives.limiter.v[0][1]);
  pj->primitives.limiter.v[1][0] =
      fmin(pi->primitives.v[1], pj->primitives.limiter.v[1][0]);
  pj->primitives.limiter.v[1][1] =
      fmax(pi->primitives.v[1], pj->primitives.limiter.v[1][1]);
  pj->primitives.limiter.v[2][0] =
      fmin(pi->primitives.v[2], pj->primitives.limiter.v[2][0]);
  pj->primitives.limiter.v[2][1] =
      fmax(pi->primitives.v[2], pj->primitives.limiter.v[2][1]);

  pj->primitives.limiter.P[0] =
      fmin(pi->primitives.P, pj->primitives.limiter.P[0]);
  pj->primitives.limiter.P[1] =
      fmax(pi->primitives.P, pj->primitives.limiter.P[1]);

  pj->primitives.limiter.maxr = fmax(r, pj->primitives.limiter.maxr);
#endif
}

/* this corresponds to task_subtype_hydro_loop1 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_hydro_loop1(float r2, float *dx, float hi, float hj,
                               struct part *pi, struct part *pj) {

  float r;
  float xi;
  float h_inv;
  float wi, wi_dx;
  int k, l;

  /* Get r and r inverse. */
  r = sqrtf(r2);

  h_inv = 1.0 / hi;
  xi = r * h_inv;
  kernel_deval(xi, &wi, &wi_dx);

  pi->density.wcount += wi;
  pi->density.wcount_dh -= xi * wi_dx;

  /* these are eqns. (1) and (2) in the summary */
  pi->geometry.volume += wi;
  for (k = 0; k < 3; k++)
    for (l = 0; l < 3; l++) pi->geometry.matrix_E[k][l] += dx[k] * dx[l] * wi;

#ifdef SPH_GRADIENTS
  /* very basic gradient estimate */
  pi->primitives.gradients.rho[0] -=
      wi_dx * dx[0] * (pi->primitives.rho - pj->primitives.rho) / r;
  pi->primitives.gradients.rho[1] -=
      wi_dx * dx[1] * (pi->primitives.rho - pj->primitives.rho) / r;
  pi->primitives.gradients.rho[2] -=
      wi_dx * dx[2] * (pi->primitives.rho - pj->primitives.rho) / r;

  pi->primitives.gradients.v[0][0] -=
      wi_dx * dx[0] * (pi->primitives.v[0] - pj->primitives.v[0]) / r;
  pi->primitives.gradients.v[0][1] -=
      wi_dx * dx[1] * (pi->primitives.v[0] - pj->primitives.v[0]) / r;
  pi->primitives.gradients.v[0][2] -=
      wi_dx * dx[2] * (pi->primitives.v[0] - pj->primitives.v[0]) / r;

  pi->primitives.gradients.v[1][0] -=
      wi_dx * dx[0] * (pi->primitives.v[1] - pj->primitives.v[1]) / r;
  pi->primitives.gradients.v[1][1] -=
      wi_dx * dx[1] * (pi->primitives.v[1] - pj->primitives.v[1]) / r;
  pi->primitives.gradients.v[1][2] -=
      wi_dx * dx[2] * (pi->primitives.v[1] - pj->primitives.v[1]) / r;

  pi->primitives.gradients.v[2][0] -=
      wi_dx * dx[0] * (pi->primitives.v[2] - pj->primitives.v[2]) / r;
  pi->primitives.gradients.v[2][1] -=
      wi_dx * dx[1] * (pi->primitives.v[2] - pj->primitives.v[2]) / r;
  pi->primitives.gradients.v[2][2] -=
      wi_dx * dx[2] * (pi->primitives.v[2] - pj->primitives.v[2]) / r;

  pi->primitives.gradients.P[0] -=
      wi_dx * dx[0] * (pi->primitives.P - pj->primitives.P) / r;
  pi->primitives.gradients.P[1] -=
      wi_dx * dx[1] * (pi->primitives.P - pj->primitives.P) / r;
  pi->primitives.gradients.P[2] -=
      wi_dx * dx[2] * (pi->primitives.P - pj->primitives.P) / r;

  /* slope limiter */
  pi->primitives.limiter.rho[0] =
      fmin(pj->primitives.rho, pi->primitives.limiter.rho[0]);
  pi->primitives.limiter.rho[1] =
      fmax(pj->primitives.rho, pi->primitives.limiter.rho[1]);

  pi->primitives.limiter.v[0][0] =
      fmin(pj->primitives.v[0], pi->primitives.limiter.v[0][0]);
  pi->primitives.limiter.v[0][1] =
      fmax(pj->primitives.v[0], pi->primitives.limiter.v[0][1]);
  pi->primitives.limiter.v[1][0] =
      fmin(pj->primitives.v[1], pi->primitives.limiter.v[1][0]);
  pi->primitives.limiter.v[1][1] =
      fmax(pj->primitives.v[1], pi->primitives.limiter.v[1][1]);
  pi->primitives.limiter.v[2][0] =
      fmin(pj->primitives.v[2], pi->primitives.limiter.v[2][0]);
  pi->primitives.limiter.v[2][1] =
      fmax(pj->primitives.v[2], pi->primitives.limiter.v[2][1]);

  pi->primitives.limiter.P[0] =
      fmin(pj->primitives.P, pi->primitives.limiter.P[0]);
  pi->primitives.limiter.P[1] =
      fmax(pj->primitives.P, pi->primitives.limiter.P[1]);

  pi->primitives.limiter.maxr = fmax(r, pi->primitives.limiter.maxr);
#endif
}

__attribute__((always_inline)) INLINE static void runner_iact_hydro_loop2(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

#ifndef SPH_GRADIENTS

  float r = sqrtf(r2);
  float xi, xj;
  float hi_inv, hj_inv;
  float wi, wj, wi_dx, wj_dx;
  int k, l;
  float Bi[3][3];
  float Bj[3][3];
  GFLOAT Wi[5], Wj[5];

  /* Initialize local variables */
  for (k = 0; k < 3; k++) {
    for (l = 0; l < 3; l++) {
      Bi[k][l] = pi->geometry.matrix_E[k][l];
      Bj[k][l] = pj->geometry.matrix_E[k][l];
    }
  }
  Wi[0] = pi->primitives.rho;
  Wi[1] = pi->primitives.v[0];
  Wi[2] = pi->primitives.v[1];
  Wi[3] = pi->primitives.v[2];
  Wi[4] = pi->primitives.P;
  Wj[0] = pj->primitives.rho;
  Wj[1] = pj->primitives.v[0];
  Wj[2] = pj->primitives.v[1];
  Wj[3] = pj->primitives.v[2];
  Wj[4] = pj->primitives.P;

  /* Compute kernel of pi. */
  hi_inv = 1.0 / hi;
  xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  /* Compute gradients for pi */
  /* there is a sign difference w.r.t. eqn. (6) because of the inverse
   * definition of dx */
  pi->primitives.gradients.rho[0] +=
      (Wi[0] - Wj[0]) * wi *
      (Bi[0][0] * dx[0] + Bi[0][1] * dx[1] + Bi[0][2] * dx[2]);
  pi->primitives.gradients.rho[1] +=
      (Wi[0] - Wj[0]) * wi *
      (Bi[1][0] * dx[0] + Bi[1][1] * dx[1] + Bi[1][2] * dx[2]);
  pi->primitives.gradients.rho[2] +=
      (Wi[0] - Wj[0]) * wi *
      (Bi[2][0] * dx[0] + Bi[2][1] * dx[1] + Bi[2][2] * dx[2]);

  pi->primitives.gradients.v[0][0] +=
      (Wi[1] - Wj[1]) * wi *
      (Bi[0][0] * dx[0] + Bi[0][1] * dx[1] + Bi[0][2] * dx[2]);
  pi->primitives.gradients.v[0][1] +=
      (Wi[1] - Wj[1]) * wi *
      (Bi[1][0] * dx[0] + Bi[1][1] * dx[1] + Bi[1][2] * dx[2]);
  pi->primitives.gradients.v[0][2] +=
      (Wi[1] - Wj[1]) * wi *
      (Bi[2][0] * dx[0] + Bi[2][1] * dx[1] + Bi[2][2] * dx[2]);
  pi->primitives.gradients.v[1][0] +=
      (Wi[2] - Wj[2]) * wi *
      (Bi[0][0] * dx[0] + Bi[0][1] * dx[1] + Bi[0][2] * dx[2]);
  pi->primitives.gradients.v[1][1] +=
      (Wi[2] - Wj[2]) * wi *
      (Bi[1][0] * dx[0] + Bi[1][1] * dx[1] + Bi[1][2] * dx[2]);
  pi->primitives.gradients.v[1][2] +=
      (Wi[2] - Wj[2]) * wi *
      (Bi[2][0] * dx[0] + Bi[2][1] * dx[1] + Bi[2][2] * dx[2]);
  pi->primitives.gradients.v[2][0] +=
      (Wi[3] - Wj[3]) * wi *
      (Bi[0][0] * dx[0] + Bi[0][1] * dx[1] + Bi[0][2] * dx[2]);
  pi->primitives.gradients.v[2][1] +=
      (Wi[3] - Wj[3]) * wi *
      (Bi[1][0] * dx[0] + Bi[1][1] * dx[1] + Bi[1][2] * dx[2]);
  pi->primitives.gradients.v[2][2] +=
      (Wi[3] - Wj[3]) * wi *
      (Bi[2][0] * dx[0] + Bi[2][1] * dx[1] + Bi[2][2] * dx[2]);

  pi->primitives.gradients.P[0] +=
      (Wi[4] - Wj[4]) * wi *
      (Bi[0][0] * dx[0] + Bi[0][1] * dx[1] + Bi[0][2] * dx[2]);
  pi->primitives.gradients.P[1] +=
      (Wi[4] - Wj[4]) * wi *
      (Bi[1][0] * dx[0] + Bi[1][1] * dx[1] + Bi[1][2] * dx[2]);
  pi->primitives.gradients.P[2] +=
      (Wi[4] - Wj[4]) * wi *
      (Bi[2][0] * dx[0] + Bi[2][1] * dx[1] + Bi[2][2] * dx[2]);

  /* basic slope limiter: collect the maximal and the minimal value for the
   * primitive variables among the ngbs */
  pi->primitives.limiter.rho[0] =
      fmin(pj->primitives.rho, pi->primitives.limiter.rho[0]);
  pi->primitives.limiter.rho[1] =
      fmax(pj->primitives.rho, pi->primitives.limiter.rho[1]);

  pi->primitives.limiter.v[0][0] =
      fmin(pj->primitives.v[0], pi->primitives.limiter.v[0][0]);
  pi->primitives.limiter.v[0][1] =
      fmax(pj->primitives.v[0], pi->primitives.limiter.v[0][1]);
  pi->primitives.limiter.v[1][0] =
      fmin(pj->primitives.v[1], pi->primitives.limiter.v[1][0]);
  pi->primitives.limiter.v[1][1] =
      fmax(pj->primitives.v[1], pi->primitives.limiter.v[1][1]);
  pi->primitives.limiter.v[2][0] =
      fmin(pj->primitives.v[2], pi->primitives.limiter.v[2][0]);
  pi->primitives.limiter.v[2][1] =
      fmax(pj->primitives.v[2], pi->primitives.limiter.v[2][1]);

  pi->primitives.limiter.P[0] =
      fmin(pj->primitives.P, pi->primitives.limiter.P[0]);
  pi->primitives.limiter.P[1] =
      fmax(pj->primitives.P, pi->primitives.limiter.P[1]);

  pi->primitives.limiter.maxr = fmax(r, pi->primitives.limiter.maxr);

  /* Compute kernel of pj. */
  hj_inv = 1.0 / hj;
  xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);

  /* Compute gradients for pj */
  /* there is no sign difference w.r.t. eqn. (6) because dx is now what we want
   * it to be */
  pj->primitives.gradients.rho[0] +=
      (Wi[0] - Wj[0]) * wj *
      (Bj[0][0] * dx[0] + Bj[0][1] * dx[1] + Bj[0][2] * dx[2]);
  pj->primitives.gradients.rho[1] +=
      (Wi[0] - Wj[0]) * wj *
      (Bj[1][0] * dx[0] + Bj[1][1] * dx[1] + Bj[1][2] * dx[2]);
  pj->primitives.gradients.rho[2] +=
      (Wi[0] - Wj[0]) * wj *
      (Bj[2][0] * dx[0] + Bj[2][1] * dx[1] + Bj[2][2] * dx[2]);

  pj->primitives.gradients.v[0][0] +=
      (Wi[1] - Wj[1]) * wj *
      (Bj[0][0] * dx[0] + Bj[0][1] * dx[1] + Bj[0][2] * dx[2]);
  pj->primitives.gradients.v[0][1] +=
      (Wi[1] - Wj[1]) * wj *
      (Bj[1][0] * dx[0] + Bj[1][1] * dx[1] + Bj[1][2] * dx[2]);
  pj->primitives.gradients.v[0][2] +=
      (Wi[1] - Wj[1]) * wj *
      (Bj[2][0] * dx[0] + Bj[2][1] * dx[1] + Bj[2][2] * dx[2]);
  pj->primitives.gradients.v[1][0] +=
      (Wi[2] - Wj[2]) * wj *
      (Bj[0][0] * dx[0] + Bj[0][1] * dx[1] + Bj[0][2] * dx[2]);
  pj->primitives.gradients.v[1][1] +=
      (Wi[2] - Wj[2]) * wj *
      (Bj[1][0] * dx[0] + Bj[1][1] * dx[1] + Bj[1][2] * dx[2]);
  pj->primitives.gradients.v[1][2] +=
      (Wi[2] - Wj[2]) * wj *
      (Bj[2][0] * dx[0] + Bj[2][1] * dx[1] + Bj[2][2] * dx[2]);
  pj->primitives.gradients.v[2][0] +=
      (Wi[3] - Wj[3]) * wj *
      (Bj[0][0] * dx[0] + Bj[0][1] * dx[1] + Bj[0][2] * dx[2]);
  pj->primitives.gradients.v[2][1] +=
      (Wi[3] - Wj[3]) * wj *
      (Bj[1][0] * dx[0] + Bj[1][1] * dx[1] + Bj[1][2] * dx[2]);
  pj->primitives.gradients.v[2][2] +=
      (Wi[3] - Wj[3]) * wj *
      (Bj[2][0] * dx[0] + Bj[2][1] * dx[1] + Bj[2][2] * dx[2]);

  pj->primitives.gradients.P[0] +=
      (Wi[4] - Wj[4]) * wj *
      (Bj[0][0] * dx[0] + Bj[0][1] * dx[1] + Bj[0][2] * dx[2]);
  pj->primitives.gradients.P[1] +=
      (Wi[4] - Wj[4]) * wj *
      (Bj[1][0] * dx[0] + Bj[1][1] * dx[1] + Bj[1][2] * dx[2]);
  pj->primitives.gradients.P[2] +=
      (Wi[4] - Wj[4]) * wj *
      (Bj[2][0] * dx[0] + Bj[2][1] * dx[1] + Bj[2][2] * dx[2]);

  /* basic slope limiter: collect the maximal and the minimal value for the
   * primitive variables among the ngbs */
  pj->primitives.limiter.rho[0] =
      fmin(pi->primitives.rho, pj->primitives.limiter.rho[0]);
  pj->primitives.limiter.rho[1] =
      fmax(pi->primitives.rho, pj->primitives.limiter.rho[1]);

  pj->primitives.limiter.v[0][0] =
      fmin(pi->primitives.v[0], pj->primitives.limiter.v[0][0]);
  pj->primitives.limiter.v[0][1] =
      fmax(pi->primitives.v[0], pj->primitives.limiter.v[0][1]);
  pj->primitives.limiter.v[1][0] =
      fmin(pi->primitives.v[1], pj->primitives.limiter.v[1][0]);
  pj->primitives.limiter.v[1][1] =
      fmax(pi->primitives.v[1], pj->primitives.limiter.v[1][1]);
  pj->primitives.limiter.v[2][0] =
      fmin(pi->primitives.v[2], pj->primitives.limiter.v[2][0]);
  pj->primitives.limiter.v[2][1] =
      fmax(pi->primitives.v[2], pj->primitives.limiter.v[2][1]);

  pj->primitives.limiter.P[0] =
      fmin(pi->primitives.P, pj->primitives.limiter.P[0]);
  pj->primitives.limiter.P[1] =
      fmax(pi->primitives.P, pj->primitives.limiter.P[1]);

  pj->primitives.limiter.maxr = fmax(r, pj->primitives.limiter.maxr);

#endif
}

__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_hydro_loop2(float r2, float *dx, float hi, float hj,
                               struct part *pi, struct part *pj) {

#ifndef SPH_GRADIENTS

  float r = sqrtf(r2);
  float xi;
  float hi_inv;
  float wi, wi_dx;
  int k, l;
  float Bi[3][3];
  GFLOAT Wi[5], Wj[5];

  /* Initialize local variables */
  for (k = 0; k < 3; k++) {
    for (l = 0; l < 3; l++) {
      Bi[k][l] = pi->geometry.matrix_E[k][l];
    }
  }
  Wi[0] = pi->primitives.rho;
  Wi[1] = pi->primitives.v[0];
  Wi[2] = pi->primitives.v[1];
  Wi[3] = pi->primitives.v[2];
  Wi[4] = pi->primitives.P;
  Wj[0] = pj->primitives.rho;
  Wj[1] = pj->primitives.v[0];
  Wj[2] = pj->primitives.v[1];
  Wj[3] = pj->primitives.v[2];
  Wj[4] = pj->primitives.P;

  /* Compute kernel of pi. */
  hi_inv = 1.0 / hi;
  xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  /* Compute gradients for pi */
  /* there is a sign difference w.r.t. eqn. (6) because of the inverse
   * definition of dx */
  pi->primitives.gradients.rho[0] +=
      (Wi[0] - Wj[0]) * wi *
      (Bi[0][0] * dx[0] + Bi[0][1] * dx[1] + Bi[0][2] * dx[2]);
  pi->primitives.gradients.rho[1] +=
      (Wi[0] - Wj[0]) * wi *
      (Bi[1][0] * dx[0] + Bi[1][1] * dx[1] + Bi[1][2] * dx[2]);
  pi->primitives.gradients.rho[2] +=
      (Wi[0] - Wj[0]) * wi *
      (Bi[2][0] * dx[0] + Bi[2][1] * dx[1] + Bi[2][2] * dx[2]);

  pi->primitives.gradients.v[0][0] +=
      (Wi[1] - Wj[1]) * wi *
      (Bi[0][0] * dx[0] + Bi[0][1] * dx[1] + Bi[0][2] * dx[2]);
  pi->primitives.gradients.v[0][1] +=
      (Wi[1] - Wj[1]) * wi *
      (Bi[1][0] * dx[0] + Bi[1][1] * dx[1] + Bi[1][2] * dx[2]);
  pi->primitives.gradients.v[0][2] +=
      (Wi[1] - Wj[1]) * wi *
      (Bi[2][0] * dx[0] + Bi[2][1] * dx[1] + Bi[2][2] * dx[2]);
  pi->primitives.gradients.v[1][0] +=
      (Wi[2] - Wj[2]) * wi *
      (Bi[0][0] * dx[0] + Bi[0][1] * dx[1] + Bi[0][2] * dx[2]);
  pi->primitives.gradients.v[1][1] +=
      (Wi[2] - Wj[2]) * wi *
      (Bi[1][0] * dx[0] + Bi[1][1] * dx[1] + Bi[1][2] * dx[2]);
  pi->primitives.gradients.v[1][2] +=
      (Wi[2] - Wj[2]) * wi *
      (Bi[2][0] * dx[0] + Bi[2][1] * dx[1] + Bi[2][2] * dx[2]);
  pi->primitives.gradients.v[2][0] +=
      (Wi[3] - Wj[3]) * wi *
      (Bi[0][0] * dx[0] + Bi[0][1] * dx[1] + Bi[0][2] * dx[2]);
  pi->primitives.gradients.v[2][1] +=
      (Wi[3] - Wj[3]) * wi *
      (Bi[1][0] * dx[0] + Bi[1][1] * dx[1] + Bi[1][2] * dx[2]);
  pi->primitives.gradients.v[2][2] +=
      (Wi[3] - Wj[3]) * wi *
      (Bi[2][0] * dx[0] + Bi[2][1] * dx[1] + Bi[2][2] * dx[2]);

  pi->primitives.gradients.P[0] +=
      (Wi[4] - Wj[4]) * wi *
      (Bi[0][0] * dx[0] + Bi[0][1] * dx[1] + Bi[0][2] * dx[2]);
  pi->primitives.gradients.P[1] +=
      (Wi[4] - Wj[4]) * wi *
      (Bi[1][0] * dx[0] + Bi[1][1] * dx[1] + Bi[1][2] * dx[2]);
  pi->primitives.gradients.P[2] +=
      (Wi[4] - Wj[4]) * wi *
      (Bi[2][0] * dx[0] + Bi[2][1] * dx[1] + Bi[2][2] * dx[2]);

  /* slope limiter */
  pi->primitives.limiter.rho[0] =
      fmin(pj->primitives.rho, pi->primitives.limiter.rho[0]);
  pi->primitives.limiter.rho[1] =
      fmax(pj->primitives.rho, pi->primitives.limiter.rho[1]);

  pi->primitives.limiter.v[0][0] =
      fmin(pj->primitives.v[0], pi->primitives.limiter.v[0][0]);
  pi->primitives.limiter.v[0][1] =
      fmax(pj->primitives.v[0], pi->primitives.limiter.v[0][1]);
  pi->primitives.limiter.v[1][0] =
      fmin(pj->primitives.v[1], pi->primitives.limiter.v[1][0]);
  pi->primitives.limiter.v[1][1] =
      fmax(pj->primitives.v[1], pi->primitives.limiter.v[1][1]);
  pi->primitives.limiter.v[2][0] =
      fmin(pj->primitives.v[2], pi->primitives.limiter.v[2][0]);
  pi->primitives.limiter.v[2][1] =
      fmax(pj->primitives.v[2], pi->primitives.limiter.v[2][1]);

  pi->primitives.limiter.P[0] =
      fmin(pj->primitives.P, pi->primitives.limiter.P[0]);
  pi->primitives.limiter.P[1] =
      fmax(pj->primitives.P, pi->primitives.limiter.P[1]);

  pi->primitives.limiter.maxr = fmax(r, pi->primitives.limiter.maxr);

#endif
}

__attribute__((always_inline)) INLINE static void runner_iact_fluxes_common(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj,
    int mode) {

  float r = sqrtf(r2);
  float xi, xj;
  float hi_inv, hi_inv2;
  float hj_inv, hj_inv2;
  float wi, wj, wi_dx, wj_dx;
  int k, l;
  float A[3];
  float Anorm;
  float Bi[3][3];
  float Bj[3][3];
  float Vi, Vj;
  float xij_i[3], xfac, xijdotdx;
  float vmax, dvdotdx;
  float vi[3], vj[3], vij[3];
  GFLOAT Wi[5], Wj[5];  //, Whalf[5];
#ifdef USE_GRADIENTS
  GFLOAT dWi[5], dWj[5];
  float xij_j[3];
#endif
  //    GFLOAT rhoe;
  //    GFLOAT flux[5][3];
  float dti, dtj, mindt;
  float n_unit[3];

  /* Initialize local variables */
  for (k = 0; k < 3; k++) {
    for (l = 0; l < 3; l++) {
      Bi[k][l] = pi->geometry.matrix_E[k][l];
      Bj[k][l] = pj->geometry.matrix_E[k][l];
    }
    vi[k] = pi->v[k]; /* particle velocities */
    vj[k] = pj->v[k];
  }
  Vi = pi->geometry.volume;
  Vj = pj->geometry.volume;
  Wi[0] = pi->primitives.rho;
  Wi[1] = pi->primitives.v[0];
  Wi[2] = pi->primitives.v[1];
  Wi[3] = pi->primitives.v[2];
  Wi[4] = pi->primitives.P;
  Wj[0] = pj->primitives.rho;
  Wj[1] = pj->primitives.v[0];
  Wj[2] = pj->primitives.v[1];
  Wj[3] = pj->primitives.v[2];
  Wj[4] = pj->primitives.P;
  dti = pi->ti_end - pi->ti_begin;  // MATTHIEU
  dtj = pj->ti_end - pj->ti_begin;

  //    if(dti > 1.e-7 || dtj > 1.e-7){
  //        message("Timestep too large: %g %g!", dti, dtj);
  //    }

  /* calculate the maximal signal velocity */
  vmax = sqrtf(const_hydro_gamma * Wi[4] / Wi[0]) +
         sqrtf(const_hydro_gamma * Wj[4] / Wj[0]);
  dvdotdx = (Wi[1] - Wj[1]) * dx[0] + (Wi[2] - Wj[2]) * dx[1] +
            (Wi[3] - Wj[3]) * dx[2];
  if (dvdotdx > 0.) {
    vmax -= dvdotdx / r;
  }
  pi->timestepvars.vmax = fmaxf(pi->timestepvars.vmax, vmax);
  if (mode == 1) {
    pj->timestepvars.vmax = fmaxf(pj->timestepvars.vmax, vmax);
  }

  /* The flux will be exchanged using the smallest time step of the two
   * particles */
  mindt = fminf(dti, dtj);

  /* Compute kernel of pi. */
  hi_inv = 1.0 / hi;
  hi_inv2 = hi_inv * hi_inv;
  xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  /* Compute kernel of pj. */
  hj_inv = 1.0 / hj;
  hj_inv2 = hj_inv * hj_inv;
  xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);

  /* Compute area */
  /* eqn. (7) */
  Anorm = 0.0f;
  for (k = 0; k < 3; k++) {
    /* we add a minus sign since dx is pi->x - pj->x */
    A[k] = -Vi * (Bi[k][0] * dx[0] + Bi[k][1] * dx[1] + Bi[k][2] * dx[2]) * wi *
               hi_inv * hi_inv2 -
           Vj * (Bj[k][0] * dx[0] + Bj[k][1] * dx[1] + Bj[k][2] * dx[2]) * wj *
               hj_inv * hj_inv2;
    Anorm += A[k] * A[k];
  }

  if (!Anorm) {
    /* if the interface has no area, nothing happens and we return */
    /* continuing results in dividing by zero and NaN's... */
    return;
  }

  /* compute the normal vector of the interface */
  Anorm = sqrtf(Anorm);
  for (k = 0; k < 3; k++) n_unit[k] = A[k] / Anorm;

#ifdef PRINT_ID
  if (pi->id == PRINT_ID || pj->id == PRINT_ID) {
    printf("pi: %g %g %g\npj: %g %g %g\nA = %g %g %g\n", pi->x[0], pi->x[1],
           pi->x[2], pj->x[0], pj->x[1], pj->x[2], A[0], A[1], A[2]);
  }
#endif

  /* Compute interface position (relative to pi, since we don't need the actual
   * position) */
  /* eqn. (8) */
  xfac = hi / (hi + hj);
  for (k = 0; k < 3; k++) xij_i[k] = -xfac * dx[k];

  /* Compute interface velocity */
  /* eqn. (9) */
  xijdotdx = xij_i[0] * dx[0] + xij_i[1] * dx[1] + xij_i[2] * dx[2];
  for (k = 0; k < 3; k++) vij[k] = vi[k] + (vi[k] - vj[k]) * xijdotdx / r2;

  /* complete calculation of position of interface */
  /* NOTE: dx is not necessarily just pi->x - pj->x but can also contain
           correction terms for periodicity. If we do the interpolation,
           we have to use xij w.r.t. the actual particle.
           => we need a separate xij for pi and pj... */
  /* tldr: we do not need the code below, but we do need the same code as above
     but then
     with i and j swapped */
  //    for ( k = 0 ; k < 3 ; k++ )
  //      xij[k] += pi->x[k];

  /* Boost the primitive variables to the frame of reference of the interface */
  /* Note that velocities are indices 1-3 in W */
  Wi[1] -= vij[0];
  Wi[2] -= vij[1];
  Wi[3] -= vij[2];
  Wj[1] -= vij[0];
  Wj[2] -= vij[1];
  Wj[3] -= vij[2];

#ifdef USE_GRADIENTS
  /* perform gradient reconstruction in space and time */
  /* space */
  /* Compute interface position (relative to pj, since we don't need the actual
   * position) */
  /* eqn. (8) */
  xfac = hj / (hi + hj);
  for (k = 0; k < 3; k++) xij_j[k] = xfac * dx[k];

  dWi[0] = pi->primitives.gradients.rho[0] * xij_i[0] +
           pi->primitives.gradients.rho[1] * xij_i[1] +
           pi->primitives.gradients.rho[2] * xij_i[2];
  dWi[1] = pi->primitives.gradients.v[0][0] * xij_i[0] +
           pi->primitives.gradients.v[0][1] * xij_i[1] +
           pi->primitives.gradients.v[0][2] * xij_i[2];
  dWi[2] = pi->primitives.gradients.v[1][0] * xij_i[0] +
           pi->primitives.gradients.v[1][1] * xij_i[1] +
           pi->primitives.gradients.v[1][2] * xij_i[2];
  dWi[3] = pi->primitives.gradients.v[2][0] * xij_i[0] +
           pi->primitives.gradients.v[2][1] * xij_i[1] +
           pi->primitives.gradients.v[2][2] * xij_i[2];
  dWi[4] = pi->primitives.gradients.P[0] * xij_i[0] +
           pi->primitives.gradients.P[1] * xij_i[1] +
           pi->primitives.gradients.P[2] * xij_i[2];

  dWj[0] = pj->primitives.gradients.rho[0] * xij_j[0] +
           pj->primitives.gradients.rho[1] * xij_j[1] +
           pj->primitives.gradients.rho[2] * xij_j[2];
  dWj[1] = pj->primitives.gradients.v[0][0] * xij_j[0] +
           pj->primitives.gradients.v[0][1] * xij_j[1] +
           pj->primitives.gradients.v[0][2] * xij_j[2];
  dWj[2] = pj->primitives.gradients.v[1][0] * xij_j[0] +
           pj->primitives.gradients.v[1][1] * xij_j[1] +
           pj->primitives.gradients.v[1][2] * xij_j[2];
  dWj[3] = pj->primitives.gradients.v[2][0] * xij_j[0] +
           pj->primitives.gradients.v[2][1] * xij_j[1] +
           pj->primitives.gradients.v[2][2] * xij_j[2];
  dWj[4] = pj->primitives.gradients.P[0] * xij_j[0] +
           pj->primitives.gradients.P[1] * xij_j[1] +
           pj->primitives.gradients.P[2] * xij_j[2];

#ifdef PER_FACE_LIMITER

  float xij_i_norm;
  GFLOAT phi_i, phi_j;
  GFLOAT delta1, delta2;
  GFLOAT phiminus, phiplus;
  GFLOAT phimin, phimax;
  GFLOAT phibar;
  /* free parameters, values from Hopkins */
  GFLOAT psi1 = 0.5, psi2 = 0.25;
  GFLOAT phi_mid0, phi_mid;

  for (k = 0; k < 10; k++) {
    if (k < 5) {
      phi_i = Wi[k];
      phi_j = Wj[k];
      phi_mid0 = Wi[k] + dWi[k];
      xij_i_norm = sqrtf(xij_i[0] * xij_i[0] + xij_i[1] * xij_i[1] +
                         xij_i[2] * xij_i[2]);
    } else {
      phi_i = Wj[k - 5];
      phi_j = Wi[k - 5];
      phi_mid0 = Wj[k - 5] + dWj[k - 5];
      xij_i_norm = sqrtf(xij_j[0] * xij_j[0] + xij_j[1] * xij_j[1] +
                         xij_j[2] * xij_j[2]);
    }

    delta1 = psi1 * fabs(phi_i - phi_j);
    delta2 = psi2 * fabs(phi_i - phi_j);

    phimin = fmin(phi_i, phi_j);
    phimax = fmax(phi_i, phi_j);

    phibar = phi_i + xij_i_norm / r * (phi_j - phi_i);

    /* if sign(phimax+delta1) == sign(phimax) */
    if ((phimax + delta1) * phimax > 0.0f) {
      phiplus = phimax + delta1;
    } else {
      phiplus = phimax / (1.0f + delta1 / fabs(phimax));
    }

    /* if sign(phimin-delta1) == sign(phimin) */
    if ((phimin - delta1) * phimin > 0.0f) {
      phiminus = phimin - delta1;
    } else {
      phiminus = phimin / (1.0f + delta1 / fabs(phimin));
    }

    if (phi_i == phi_j) {
      phi_mid = phi_i;
    } else {
      if (phi_i < phi_j) {
        phi_mid = fmax(phiminus, fmin(phibar + delta2, phi_mid0));
      } else {
        phi_mid = fmin(phiplus, fmax(phibar - delta2, phi_mid0));
      }
    }

    if (k < 5) {
      dWi[k] = phi_mid - phi_i;
    } else {
      dWj[k - 5] = phi_mid - phi_i;
    }
  }

#endif

  //    printf("dWL: %g %g %g %g %g\n", dWi[0], dWi[1], dWi[2], dWi[3], dWi[4]);
  //    printf("dWR: %g %g %g %g %g\n", dWj[0], dWj[1], dWj[2], dWj[3], dWj[4]);

  /* time */
  dWi[0] -= 0.5 * mindt * (Wi[1] * pi->primitives.gradients.rho[0] +
                           Wi[2] * pi->primitives.gradients.rho[1] +
                           Wi[3] * pi->primitives.gradients.rho[2] +
                           Wi[0] * (pi->primitives.gradients.v[0][0] +
                                    pi->primitives.gradients.v[1][1] +
                                    pi->primitives.gradients.v[2][2]));
  dWi[1] -= 0.5 * mindt * (Wi[1] * pi->primitives.gradients.v[0][0] +
                           Wi[2] * pi->primitives.gradients.v[0][1] +
                           Wi[3] * pi->primitives.gradients.v[0][2] +
                           pi->primitives.gradients.P[0] / Wi[0]);
  dWi[2] -= 0.5 * mindt * (Wi[1] * pi->primitives.gradients.v[1][0] +
                           Wi[2] * pi->primitives.gradients.v[1][1] +
                           Wi[3] * pi->primitives.gradients.v[1][2] +
                           pi->primitives.gradients.P[1] / Wi[0]);
  dWi[3] -= 0.5 * mindt * (Wi[1] * pi->primitives.gradients.v[2][0] +
                           Wi[2] * pi->primitives.gradients.v[2][1] +
                           Wi[3] * pi->primitives.gradients.v[2][2] +
                           pi->primitives.gradients.P[2] / Wi[0]);
  dWi[4] -= 0.5 * mindt *
            (Wi[1] * pi->primitives.gradients.P[0] +
             Wi[2] * pi->primitives.gradients.P[1] +
             Wi[3] * pi->primitives.gradients.P[2] +
             const_hydro_gamma * Wi[4] * (pi->primitives.gradients.v[0][0] +
                                          pi->primitives.gradients.v[1][1] +
                                          pi->primitives.gradients.v[2][2]));

  dWj[0] -= 0.5 * mindt * (Wj[1] * pj->primitives.gradients.rho[0] +
                           Wj[2] * pj->primitives.gradients.rho[1] +
                           Wj[3] * pj->primitives.gradients.rho[2] +
                           Wj[0] * (pj->primitives.gradients.v[0][0] +
                                    pj->primitives.gradients.v[1][1] +
                                    pj->primitives.gradients.v[2][2]));
  dWj[1] -= 0.5 * mindt * (Wj[1] * pj->primitives.gradients.v[0][0] +
                           Wj[2] * pj->primitives.gradients.v[0][1] +
                           Wj[3] * pj->primitives.gradients.v[0][2] +
                           pj->primitives.gradients.P[0] / Wj[0]);
  dWj[2] -= 0.5 * mindt * (Wj[1] * pj->primitives.gradients.v[1][0] +
                           Wj[2] * pj->primitives.gradients.v[1][1] +
                           Wj[3] * pj->primitives.gradients.v[1][2] +
                           pj->primitives.gradients.P[1] / Wj[0]);
  dWj[3] -= 0.5 * mindt * (Wj[1] * pj->primitives.gradients.v[2][0] +
                           Wj[2] * pj->primitives.gradients.v[2][1] +
                           Wj[3] * pj->primitives.gradients.v[2][2] +
                           pj->primitives.gradients.P[2] / Wj[0]);
  dWj[4] -= 0.5 * mindt *
            (Wj[1] * pj->primitives.gradients.P[0] +
             Wj[2] * pj->primitives.gradients.P[1] +
             Wj[3] * pj->primitives.gradients.P[2] +
             const_hydro_gamma * Wj[4] * (pj->primitives.gradients.v[0][0] +
                                          pj->primitives.gradients.v[1][1] +
                                          pj->primitives.gradients.v[2][2]));

  //    printf("WL: %g %g %g %g %g\n", Wi[0], Wi[1], Wi[2], Wi[3], Wi[4]);
  //    printf("WR: %g %g %g %g %g\n", Wj[0], Wj[1], Wj[2], Wj[3], Wj[4]);

  //    printf("dWL: %g %g %g %g %g\n", dWi[0], dWi[1], dWi[2], dWi[3], dWi[4]);
  //    printf("dWR: %g %g %g %g %g\n", dWj[0], dWj[1], dWj[2], dWj[3], dWj[4]);

  Wi[0] += dWi[0];
  Wi[1] += dWi[1];
  Wi[2] += dWi[2];
  Wi[3] += dWi[3];
  Wi[4] += dWi[4];

  Wj[0] += dWj[0];
  Wj[1] += dWj[1];
  Wj[2] += dWj[2];
  Wj[3] += dWj[3];
  Wj[4] += dWj[4];
#endif

  /* apply slope limiter interface by interface */
  /* ... to be done ... */

  /* we don't need to rotate, we can use the unit vector in the Riemann problem
   * itself (see GIZMO) */

  if (Wi[0] < 0.0f || Wj[0] < 0.0f || Wi[4] < 0.0f || Wj[4] < 0.0f) {
    printf("mindt: %g\n", mindt);
    printf("WL: %g %g %g %g %g\n", pi->primitives.rho, pi->primitives.v[0],
           pi->primitives.v[1], pi->primitives.v[2], pi->primitives.P);
    printf("dWL: %g %g %g %g %g\n", dWi[0], dWi[1], dWi[2], dWi[3], dWi[4]);
    printf("gradWL[0]: %g %g %g\n", pi->primitives.gradients.rho[0],
           pi->primitives.gradients.rho[1], pi->primitives.gradients.rho[2]);
    printf("gradWL[1]: %g %g %g\n", pi->primitives.gradients.v[0][0],
           pi->primitives.gradients.v[0][1], pi->primitives.gradients.v[0][2]);
    printf("gradWL[2]: %g %g %g\n", pi->primitives.gradients.v[1][0],
           pi->primitives.gradients.v[1][1], pi->primitives.gradients.v[1][2]);
    printf("gradWL[3]: %g %g %g\n", pi->primitives.gradients.v[2][0],
           pi->primitives.gradients.v[2][1], pi->primitives.gradients.v[2][2]);
    printf("gradWL[4]: %g %g %g\n", pi->primitives.gradients.P[0],
           pi->primitives.gradients.P[1], pi->primitives.gradients.P[2]);
    printf("WL': %g %g %g %g %g\n", Wi[0], Wi[1], Wi[2], Wi[3], Wi[4]);
    printf("WR: %g %g %g %g %g\n", pj->primitives.rho, pj->primitives.v[0],
           pj->primitives.v[1], pj->primitives.v[2], pj->primitives.P);
    printf("dWR: %g %g %g %g %g\n", dWj[0], dWj[1], dWj[2], dWj[3], dWj[4]);
    printf("gradWR[0]: %g %g %g\n", pj->primitives.gradients.rho[0],
           pj->primitives.gradients.rho[1], pj->primitives.gradients.rho[2]);
    printf("gradWR[1]: %g %g %g\n", pj->primitives.gradients.v[0][0],
           pj->primitives.gradients.v[0][1], pj->primitives.gradients.v[0][2]);
    printf("gradWR[2]: %g %g %g\n", pj->primitives.gradients.v[1][0],
           pj->primitives.gradients.v[1][1], pj->primitives.gradients.v[1][2]);
    printf("gradWR[3]: %g %g %g\n", pj->primitives.gradients.v[2][0],
           pj->primitives.gradients.v[2][1], pj->primitives.gradients.v[2][2]);
    printf("gradWR[4]: %g %g %g\n", pj->primitives.gradients.P[0],
           pj->primitives.gradients.P[1], pj->primitives.gradients.P[2]);
    printf("WR': %g %g %g %g %g\n", Wj[0], Wj[1], Wj[2], Wj[3], Wj[4]);
    error("Negative density or pressure!\n");
  }

  GFLOAT totflux[5];
  riemann_solve_for_flux(Wi, Wj, n_unit, vij, totflux);

  /* Update conserved variables */
  /* eqn. (16) */
  pi->conserved.mass -= mindt * Anorm * totflux[0];
  pi->conserved.momentum[0] -= mindt * Anorm * totflux[1];
  pi->conserved.momentum[1] -= mindt * Anorm * totflux[2];
  pi->conserved.momentum[2] -= mindt * Anorm * totflux[3];
  pi->conserved.energy -= mindt * Anorm * totflux[4];

#ifdef THERMAL_ENERGY
  float ekin = 0.5 * (pi->primitives.v[0] * pi->primitives.v[0] +
                      pi->primitives.v[1] * pi->primitives.v[1] +
                      pi->primitives.v[2] * pi->primitives.v[2]);
  pi->conserved.energy += mindt * Anorm * totflux[1] * pi->primitives.v[0];
  pi->conserved.energy += mindt * Anorm * totflux[2] * pi->primitives.v[1];
  pi->conserved.energy += mindt * Anorm * totflux[3] * pi->primitives.v[2];
  pi->conserved.energy -= mindt * Anorm * totflux[0] * ekin;
#endif

  /* the non symmetric version is never called when using mindt, whether this
   * piece of code
   * should always be executed or only in the symmetric case is currently
   * unclear */
  if (mode == 1) {
    pj->conserved.mass += mindt * Anorm * totflux[0];
    pj->conserved.momentum[0] += mindt * Anorm * totflux[1];
    pj->conserved.momentum[1] += mindt * Anorm * totflux[2];
    pj->conserved.momentum[2] += mindt * Anorm * totflux[3];
    pj->conserved.energy += mindt * Anorm * totflux[4];

#ifdef THERMAL_ENERGY
    ekin = 0.5 * (pj->primitives.v[0] * pj->primitives.v[0] +
                  pj->primitives.v[1] * pj->primitives.v[1] +
                  pj->primitives.v[2] * pj->primitives.v[2]);
    pj->conserved.energy -= mindt * Anorm * totflux[1] * pj->primitives.v[0];
    pj->conserved.energy -= mindt * Anorm * totflux[2] * pj->primitives.v[1];
    pj->conserved.energy -= mindt * Anorm * totflux[3] * pj->primitives.v[2];
    pj->conserved.energy += mindt * Anorm * totflux[0] * ekin;
#endif
  }
}

/* this corresponds to task_subtype_fluxes */
__attribute__((always_inline)) INLINE static void runner_iact_hydro_loop3(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  runner_iact_fluxes_common(r2, dx, hi, hj, pi, pj, 1);
}

/* this corresponds to task_subtype_fluxes */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_hydro_loop3(float r2, float *dx, float hi, float hj,
                               struct part *pi, struct part *pj) {

  runner_iact_fluxes_common(r2, dx, hi, hj, pi, pj, 0);
}
