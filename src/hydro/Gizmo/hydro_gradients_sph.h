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

  p->primitives.limiter.rho[0] = FLT_MAX;
  p->primitives.limiter.rho[1] = -FLT_MAX;
  p->primitives.limiter.v[0][0] = FLT_MAX;
  p->primitives.limiter.v[0][1] = -FLT_MAX;
  p->primitives.limiter.v[1][0] = FLT_MAX;
  p->primitives.limiter.v[1][1] = -FLT_MAX;
  p->primitives.limiter.v[2][0] = FLT_MAX;
  p->primitives.limiter.v[2][1] = -FLT_MAX;
  p->primitives.limiter.P[0] = FLT_MAX;
  p->primitives.limiter.P[1] = -FLT_MAX;

  p->primitives.limiter.maxr = -FLT_MAX;
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
  }
}

/**
 * @brief Calculations done before the force loop
 */
__attribute__((always_inline)) INLINE static void
hydro_gradients_prepare_force_loop(struct part *p, float ih2, float volume) {

  float gradrho[3], gradv[3][3], gradP[3];
  float gradtrue, gradmax, gradmin, alpha;

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

  /* slope limiter */
  gradrho[0] = p->primitives.gradients.rho[0];
  gradrho[1] = p->primitives.gradients.rho[1];
  gradrho[2] = p->primitives.gradients.rho[2];

  gradv[0][0] = p->primitives.gradients.v[0][0];
  gradv[0][1] = p->primitives.gradients.v[0][1];
  gradv[0][2] = p->primitives.gradients.v[0][2];

  gradv[1][0] = p->primitives.gradients.v[1][0];
  gradv[1][1] = p->primitives.gradients.v[1][1];
  gradv[1][2] = p->primitives.gradients.v[1][2];

  gradv[2][0] = p->primitives.gradients.v[2][0];
  gradv[2][1] = p->primitives.gradients.v[2][1];
  gradv[2][2] = p->primitives.gradients.v[2][2];

  gradP[0] = p->primitives.gradients.P[0];
  gradP[1] = p->primitives.gradients.P[1];
  gradP[2] = p->primitives.gradients.P[2];

  gradtrue = sqrtf(gradrho[0] * gradrho[0] + gradrho[1] * gradrho[1] +
                   gradrho[2] * gradrho[2]); /* gradtrue might be zero. In this
    case, there is no gradient and we don't
    need to slope limit anything... */
  if (gradtrue) {
    gradtrue *= p->primitives.limiter.maxr;
    gradmax = p->primitives.limiter.rho[1] - p->primitives.rho;
    gradmin = p->primitives.rho - p->primitives.limiter.rho[0];
    alpha = fmin(1.0f, fmin(gradmax / gradtrue, gradmin / gradtrue));
    p->primitives.gradients.rho[0] *= alpha;
    p->primitives.gradients.rho[1] *= alpha;
    p->primitives.gradients.rho[2] *= alpha;
  }

  gradtrue = sqrtf(gradv[0][0] * gradv[0][0] + gradv[0][1] * gradv[0][1] +
                   gradv[0][2] * gradv[0][2]);
  if (gradtrue) {
    gradtrue *= p->primitives.limiter.maxr;
    gradmax = p->primitives.limiter.v[0][1] - p->primitives.v[0];
    gradmin = p->primitives.v[0] - p->primitives.limiter.v[0][0];
    alpha = fmin(1.0f, fmin(gradmax / gradtrue, gradmin / gradtrue));
    p->primitives.gradients.v[0][0] *= alpha;
    p->primitives.gradients.v[0][1] *= alpha;
    p->primitives.gradients.v[0][2] *= alpha;
  }

  gradtrue = sqrtf(gradv[1][0] * gradv[1][0] + gradv[1][1] * gradv[1][1] +
                   gradv[1][2] * gradv[1][2]);
  if (gradtrue) {
    gradtrue *= p->primitives.limiter.maxr;
    gradmax = p->primitives.limiter.v[1][1] - p->primitives.v[1];
    gradmin = p->primitives.v[1] - p->primitives.limiter.v[1][0];
    alpha = fmin(1.0f, fmin(gradmax / gradtrue, gradmin / gradtrue));
    p->primitives.gradients.v[1][0] *= alpha;
    p->primitives.gradients.v[1][1] *= alpha;
    p->primitives.gradients.v[1][2] *= alpha;
  }

  gradtrue = sqrtf(gradv[2][0] * gradv[2][0] + gradv[2][1] * gradv[2][1] +
                   gradv[2][2] * gradv[2][2]);
  if (gradtrue) {
    gradtrue *= p->primitives.limiter.maxr;
    gradmax = p->primitives.limiter.v[2][1] - p->primitives.v[2];
    gradmin = p->primitives.v[2] - p->primitives.limiter.v[2][0];
    alpha = fmin(1.0f, fmin(gradmax / gradtrue, gradmin / gradtrue));
    p->primitives.gradients.v[2][0] *= alpha;
    p->primitives.gradients.v[2][1] *= alpha;
    p->primitives.gradients.v[2][2] *= alpha;
  }

  gradtrue =
      sqrtf(gradP[0] * gradP[0] + gradP[1] * gradP[1] + gradP[2] * gradP[2]);
  if (gradtrue) {
    gradtrue *= p->primitives.limiter.maxr;
    gradmax = p->primitives.limiter.P[1] - p->primitives.P;
    gradmin = p->primitives.P - p->primitives.limiter.P[0];
    alpha = fmin(1.0f, fmin(gradmax / gradtrue, gradmin / gradtrue));
    p->primitives.gradients.P[0] *= alpha;
    p->primitives.gradients.P[1] *= alpha;
    p->primitives.gradients.P[2] *= alpha;
  }
}

/**
 * @brief Gradient calculations done during the gradient loop
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_gradient_loop(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj,
    int mode) {}
