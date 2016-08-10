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
 * @brief Initialize variables before the density loop
 */
__attribute__((always_inline)) INLINE static void
hydro_gradients_init_density_loop(struct part *p) {

  /* use the old volumes to estimate new primitive variables to be used for the
     gradient calculation */
  if (p->conserved.mass) {
    p->primitives.rho = p->conserved.mass / p->geometry.volume;
    p->primitives.v[0] = p->conserved.momentum[0] / p->conserved.mass;
    p->primitives.v[1] = p->conserved.momentum[1] / p->conserved.mass;
    p->primitives.v[2] = p->conserved.momentum[2] / p->conserved.mass;
    p->primitives.P =
        hydro_gamma_minus_one *
        (p->conserved.energy -
         0.5 * (p->conserved.momentum[0] * p->conserved.momentum[0] +
                p->conserved.momentum[1] * p->conserved.momentum[1] +
                p->conserved.momentum[2] * p->conserved.momentum[2]) /
             p->conserved.mass) /
        p->geometry.volume;
  }

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
 * @brief Gradient calculations done during the density loop
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_density_loop(
    struct part *pi, struct part *pj, float wi_dx, float wj_dx, float *dx,
    float r, int mode) {

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

  hydro_slope_limit_cell_collect(pi, pj, r);

  if (mode == 1) {
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

    hydro_slope_limit_cell_collect(pj, pi, r);
  }
}

/**
 * @brief Calculations done before the force loop
 */
__attribute__((always_inline)) INLINE static void
hydro_gradients_prepare_force_loop(struct part *p, float ih2, float volume) {

  /* finalize gradients by multiplying with volume */
  p->primitives.gradients.rho[0] *= ih2 * ih2 * volume;
  p->primitives.gradients.rho[1] *= ih2 * ih2 * volume;
  p->primitives.gradients.rho[2] *= ih2 * ih2 * volume;

  p->primitives.gradients.v[0][0] *= ih2 * ih2 * volume;
  p->primitives.gradients.v[0][1] *= ih2 * ih2 * volume;
  p->primitives.gradients.v[0][2] *= ih2 * ih2 * volume;

  p->primitives.gradients.v[1][0] *= ih2 * ih2 * volume;
  p->primitives.gradients.v[1][1] *= ih2 * ih2 * volume;
  p->primitives.gradients.v[1][2] *= ih2 * ih2 * volume;

  p->primitives.gradients.v[2][0] *= ih2 * ih2 * volume;
  p->primitives.gradients.v[2][1] *= ih2 * ih2 * volume;
  p->primitives.gradients.v[2][2] *= ih2 * ih2 * volume;

  p->primitives.gradients.P[0] *= ih2 * ih2 * volume;
  p->primitives.gradients.P[1] *= ih2 * ih2 * volume;
  p->primitives.gradients.P[2] *= ih2 * ih2 * volume;

  hydro_slope_limit_cell(p);
}

/**
 * @brief Gradient calculations done during the gradient loop
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_gradient_loop(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj,
    int mode) {}
