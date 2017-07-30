
/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

#ifndef SWIFT_HYDRO_GRADIENTS_H
#define SWIFT_HYDRO_GRADIENTS_H

#include "hydro_slope_limiters.h"
#include "hydro_unphysical.h"
#include "riemann.h"

#if defined(GRADIENTS_SPH)

#define HYDRO_GRADIENT_IMPLEMENTATION "SPH gradients (Price 2012)"
#include "hydro_gradients_sph.h"

#elif defined(GRADIENTS_GIZMO)

#define HYDRO_GRADIENT_IMPLEMENTATION "GIZMO gradients (Hopkins 2015)"
#include "hydro_gradients_gizmo.h"

#else

/* No gradients. Perfectly acceptable, but we have to provide empty functions */
#define HYDRO_GRADIENT_IMPLEMENTATION "No gradients (first order scheme)"

/**
 * @brief Initialize gradient variables
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_init(
    struct part* p) {}

/**
 * @brief Gradient calculations done during the neighbour loop
 *
 * @param r2 Squared distance between the two particles.
 * @param dx Distance vector (pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_collect(
    float r2, float* dx, float hi, float hj, struct part* pi, struct part* pj) {
}

/**
 * @brief Gradient calculations done during the neighbour loop: non-symmetric
 * version
 *
 * @param r2 Squared distance between the two particles.
 * @param dx Distance vector (pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void
hydro_gradients_nonsym_collect(float r2, float* dx, float hi, float hj,
                               struct part* pi, struct part* pj) {}

/**
 * @brief Finalize the gradient variables after all data have been collected
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_finalize(
    struct part* p) {}

#endif

/**
 * @brief Gradients reconstruction. Is the same for all gradient types (although
 * gradients_none does nothing, since all gradients are zero -- are they?).
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_predict(
    struct part* pi, struct part* pj, float hi, float hj, float* dx, float r,
    float* xij_i, float* Wi, float* Wj, float mindt) {

  float dWi[5], dWj[5];
  float xij_j[3];
  int k;
  float xfac;

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

  hydro_slope_limit_face(Wi, Wj, dWi, dWj, xij_i, xij_j, r);

  /* time */
  if (Wi[0] > 0.0f) {
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
               hydro_gamma * Wi[4] * (pi->primitives.gradients.v[0][0] +
                                      pi->primitives.gradients.v[1][1] +
                                      pi->primitives.gradients.v[2][2]));
  }

  if (Wj[0] > 0.0f) {
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
               hydro_gamma * Wj[4] * (pj->primitives.gradients.v[0][0] +
                                      pj->primitives.gradients.v[1][1] +
                                      pj->primitives.gradients.v[2][2]));
  }

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

  gizmo_check_physical_quantity("density", Wi[0]);
  gizmo_check_physical_quantity("pressure", Wi[4]);
  gizmo_check_physical_quantity("density", Wj[0]);
  gizmo_check_physical_quantity("pressure", Wj[4]);
}

#endif  // SWIFT_HYDRO_GRADIENTS_H
