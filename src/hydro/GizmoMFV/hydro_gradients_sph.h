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

/**
 * @brief Initialize gradient variables
 *
 * @param p Particle.
 */
#ifndef SWIFT_GIZMO_MFV_HYDRO_SPH_GRADIENTS_H
#define SWIFT_GIZMO_MFV_HYDRO_SPH_GRADIENTS_H

__attribute__((always_inline)) INLINE static void hydro_gradients_init(
    struct part *p) {

  p->primitives.gradients.rho[0] = 0.0f;
  p->primitives.gradients.rho[1] = 0.0f;
  p->primitives.gradients.rho[2] = 0.0f;

  p->primitives.gradients.v[0][0] = 0.0f;
  p->primitives.gradients.v[0][1] = 0.0f;
  p->primitives.gradients.v[0][2] = 0.0f;

  p->primitives.gradients.v[1][0] = 0.0f;
  p->primitives.gradients.v[1][1] = 0.0f;
  p->primitives.gradients.v[1][2] = 0.0f;
  p->primitives.gradients.v[2][0] = 0.0f;
  p->primitives.gradients.v[2][1] = 0.0f;
  p->primitives.gradients.v[2][2] = 0.0f;

  p->primitives.gradients.P[0] = 0.0f;
  p->primitives.gradients.P[1] = 0.0f;
  p->primitives.gradients.P[2] = 0.0f;

  hydro_slope_limit_cell_init(p);
}

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
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj) {

  const float r_inv = 1.f / sqrtf(r2);
  const float r = r2 * r_inv;

  float wi, wi_dx;
  const float hi_inv = 1.0f / hi;
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  /* very basic gradient estimate */
  pi->primitives.gradients.rho[0] -=
      wi_dx * dx[0] * (pi->primitives.rho - pj->primitives.rho) * r_inv;
  pi->primitives.gradients.rho[1] -=
      wi_dx * dx[1] * (pi->primitives.rho - pj->primitives.rho) * r_inv;
  pi->primitives.gradients.rho[2] -=
      wi_dx * dx[2] * (pi->primitives.rho - pj->primitives.rho) * r_inv;

  pi->primitives.gradients.v[0][0] -=
      wi_dx * dx[0] * (pi->primitives.v[0] - pj->primitives.v[0]) * r_inv;
  pi->primitives.gradients.v[0][1] -=
      wi_dx * dx[1] * (pi->primitives.v[0] - pj->primitives.v[0]) * r_inv;
  pi->primitives.gradients.v[0][2] -=
      wi_dx * dx[2] * (pi->primitives.v[0] - pj->primitives.v[0]) * r_inv;

  pi->primitives.gradients.v[1][0] -=
      wi_dx * dx[0] * (pi->primitives.v[1] - pj->primitives.v[1]) * r_inv;
  pi->primitives.gradients.v[1][1] -=
      wi_dx * dx[1] * (pi->primitives.v[1] - pj->primitives.v[1]) * r_inv;
  pi->primitives.gradients.v[1][2] -=
      wi_dx * dx[2] * (pi->primitives.v[1] - pj->primitives.v[1]) * r_inv;

  pi->primitives.gradients.v[2][0] -=
      wi_dx * dx[0] * (pi->primitives.v[2] - pj->primitives.v[2]) * r_inv;
  pi->primitives.gradients.v[2][1] -=
      wi_dx * dx[1] * (pi->primitives.v[2] - pj->primitives.v[2]) * r_inv;
  pi->primitives.gradients.v[2][2] -=
      wi_dx * dx[2] * (pi->primitives.v[2] - pj->primitives.v[2]) * r_inv;

  pi->primitives.gradients.P[0] -=
      wi_dx * dx[0] * (pi->primitives.P - pj->primitives.P) * r_inv;
  pi->primitives.gradients.P[1] -=
      wi_dx * dx[1] * (pi->primitives.P - pj->primitives.P) * r_inv;
  pi->primitives.gradients.P[2] -=
      wi_dx * dx[2] * (pi->primitives.P - pj->primitives.P) * r_inv;

  hydro_slope_limit_cell_collect(pi, pj, r);

  float wj, wj_dx;
  const float hj_inv = 1.0f / hj;
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);

  /* signs are the same as before, since we swap i and j twice */
  pj->primitives.gradients.rho[0] -=
      wj_dx * dx[0] * (pi->primitives.rho - pj->primitives.rho) * r_inv;
  pj->primitives.gradients.rho[1] -=
      wj_dx * dx[1] * (pi->primitives.rho - pj->primitives.rho) * r_inv;
  pj->primitives.gradients.rho[2] -=
      wj_dx * dx[2] * (pi->primitives.rho - pj->primitives.rho) * r_inv;

  pj->primitives.gradients.v[0][0] -=
      wj_dx * dx[0] * (pi->primitives.v[0] - pj->primitives.v[0]) * r_inv;
  pj->primitives.gradients.v[0][1] -=
      wj_dx * dx[1] * (pi->primitives.v[0] - pj->primitives.v[0]) * r_inv;
  pj->primitives.gradients.v[0][2] -=
      wj_dx * dx[2] * (pi->primitives.v[0] - pj->primitives.v[0]) * r_inv;

  pj->primitives.gradients.v[1][0] -=
      wj_dx * dx[0] * (pi->primitives.v[1] - pj->primitives.v[1]) * r_inv;
  pj->primitives.gradients.v[1][1] -=
      wj_dx * dx[1] * (pi->primitives.v[1] - pj->primitives.v[1]) * r_inv;
  pj->primitives.gradients.v[1][2] -=
      wj_dx * dx[2] * (pi->primitives.v[1] - pj->primitives.v[1]) * r_inv;
  pj->primitives.gradients.v[2][0] -=
      wj_dx * dx[0] * (pi->primitives.v[2] - pj->primitives.v[2]) * r_inv;
  pj->primitives.gradients.v[2][1] -=
      wj_dx * dx[1] * (pi->primitives.v[2] - pj->primitives.v[2]) * r_inv;
  pj->primitives.gradients.v[2][2] -=
      wj_dx * dx[2] * (pi->primitives.v[2] - pj->primitives.v[2]) * r_inv;

  pj->primitives.gradients.P[0] -=
      wj_dx * dx[0] * (pi->primitives.P - pj->primitives.P) * r_inv;
  pj->primitives.gradients.P[1] -=
      wj_dx * dx[1] * (pi->primitives.P - pj->primitives.P) * r_inv;
  pj->primitives.gradients.P[2] -=
      wj_dx * dx[2] * (pi->primitives.P - pj->primitives.P) * r_inv;

  hydro_slope_limit_cell_collect(pj, pi, r);
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
hydro_gradients_nonsym_collect(float r2, const float *dx, float hi, float hj,
                               struct part *restrict pi,
                               struct part *restrict pj) {

  const float r_inv = 1.f / sqrtf(r2);
  const float r = r2 * r_inv;

  float wi, wi_dx;
  const float hi_inv = 1.0f / hi;
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  /* very basic gradient estimate */
  pi->primitives.gradients.rho[0] -=
      wi_dx * dx[0] * (pi->primitives.rho - pj->primitives.rho) * r_inv;
  pi->primitives.gradients.rho[1] -=
      wi_dx * dx[1] * (pi->primitives.rho - pj->primitives.rho) * r_inv;
  pi->primitives.gradients.rho[2] -=
      wi_dx * dx[2] * (pi->primitives.rho - pj->primitives.rho) * r_inv;

  pi->primitives.gradients.v[0][0] -=
      wi_dx * dx[0] * (pi->primitives.v[0] - pj->primitives.v[0]) * r_inv;
  pi->primitives.gradients.v[0][1] -=
      wi_dx * dx[1] * (pi->primitives.v[0] - pj->primitives.v[0]) * r_inv;
  pi->primitives.gradients.v[0][2] -=
      wi_dx * dx[2] * (pi->primitives.v[0] - pj->primitives.v[0]) * r_inv;

  pi->primitives.gradients.v[1][0] -=
      wi_dx * dx[0] * (pi->primitives.v[1] - pj->primitives.v[1]) * r_inv;
  pi->primitives.gradients.v[1][1] -=
      wi_dx * dx[1] * (pi->primitives.v[1] - pj->primitives.v[1]) * r_inv;
  pi->primitives.gradients.v[1][2] -=
      wi_dx * dx[2] * (pi->primitives.v[1] - pj->primitives.v[1]) * r_inv;

  pi->primitives.gradients.v[2][0] -=
      wi_dx * dx[0] * (pi->primitives.v[2] - pj->primitives.v[2]) * r_inv;
  pi->primitives.gradients.v[2][1] -=
      wi_dx * dx[1] * (pi->primitives.v[2] - pj->primitives.v[2]) * r_inv;
  pi->primitives.gradients.v[2][2] -=
      wi_dx * dx[2] * (pi->primitives.v[2] - pj->primitives.v[2]) * r_inv;

  pi->primitives.gradients.P[0] -=
      wi_dx * dx[0] * (pi->primitives.P - pj->primitives.P) * r_inv;
  pi->primitives.gradients.P[1] -=
      wi_dx * dx[1] * (pi->primitives.P - pj->primitives.P) * r_inv;
  pi->primitives.gradients.P[2] -=
      wi_dx * dx[2] * (pi->primitives.P - pj->primitives.P) * r_inv;

  hydro_slope_limit_cell_collect(pi, pj, r);
}

/**
 * @brief Finalize the gradient variables after all data have been collected
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_finalize(
    struct part *p) {

  const float h = p->h;
  const float ih = 1.0f / h;
  const float ihdimp1 = pow_dimension_plus_one(ih);
  const float volume = p->geometry.volume;

  /* finalize gradients by multiplying with volume */
  p->primitives.gradients.rho[0] *= ihdimp1 * volume;
  p->primitives.gradients.rho[1] *= ihdimp1 * volume;
  p->primitives.gradients.rho[2] *= ihdimp1 * volume;

  p->primitives.gradients.v[0][0] *= ihdimp1 * volume;
  p->primitives.gradients.v[0][1] *= ihdimp1 * volume;
  p->primitives.gradients.v[0][2] *= ihdimp1 * volume;

  p->primitives.gradients.v[1][0] *= ihdimp1 * volume;
  p->primitives.gradients.v[1][1] *= ihdimp1 * volume;
  p->primitives.gradients.v[1][2] *= ihdimp1 * volume;

  p->primitives.gradients.v[2][0] *= ihdimp1 * volume;
  p->primitives.gradients.v[2][1] *= ihdimp1 * volume;
  p->primitives.gradients.v[2][2] *= ihdimp1 * volume;

  p->primitives.gradients.P[0] *= ihdimp1 * volume;
  p->primitives.gradients.P[1] *= ihdimp1 * volume;
  p->primitives.gradients.P[2] *= ihdimp1 * volume;

  hydro_slope_limit_cell(p);
}

#endif /* SWIFT_GIZMO_MFV_HYDRO_SPH_GRADIENTS_H */
