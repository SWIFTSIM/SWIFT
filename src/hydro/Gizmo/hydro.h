/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
 * @brief Computes the hydro time-step of a given particle
 *
 * @param p Pointer to the particle data
 * @param xp Pointer to the extended particle data
 *
 */
__attribute__((always_inline)) INLINE static float hydro_compute_timestep(
    struct part* p, struct xpart* xp) {

  return const_cfl * p->h / fabs(p->timestepvars.vmax);
}

/**
 * @brief Initialises the particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_first_init_part(
    struct part* p, struct xpart* xp) {}

/**
 * @brief Prepares a particle for the volume calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_init_part(
    struct part* p) {

#ifdef SPH_GRADIENTS
  /* use the old volumes to estimate new primitive variables to be used for the
     gradient calculation */
  if (p->conserved.mass) {
    p->primitives.rho = p->conserved.mass / p->geometry.volume;
    p->primitives.v[0] = p->conserved.momentum[0] / p->conserved.mass;
    p->primitives.v[1] = p->conserved.momentum[1] / p->conserved.mass;
    p->primitives.v[2] = p->conserved.momentum[2] / p->conserved.mass;
    p->primitives.P =
        (const_hydro_gamma - 1.) *
        (p->conserved.energy -
         0.5 * (p->conserved.momentum[0] * p->conserved.momentum[0] +
                p->conserved.momentum[1] * p->conserved.momentum[1] +
                p->conserved.momentum[2] * p->conserved.momentum[2]) /
             p->conserved.mass) /
        p->geometry.volume;
  }
#endif

  p->density.wcount = 0.0f;
  p->density.wcount_dh = 0.0f;
  p->geometry.volume = 0.0f;
  p->geometry.matrix_E[0][0] = 0.0f;
  p->geometry.matrix_E[0][1] = 0.0f;
  p->geometry.matrix_E[0][2] = 0.0f;
  p->geometry.matrix_E[1][0] = 0.0f;
  p->geometry.matrix_E[1][1] = 0.0f;
  p->geometry.matrix_E[1][2] = 0.0f;
  p->geometry.matrix_E[2][0] = 0.0f;
  p->geometry.matrix_E[2][1] = 0.0f;
  p->geometry.matrix_E[2][2] = 0.0f;

#ifdef SPH_GRADIENTS
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
#endif
}

/**
 * @brief Finishes the density calculation.
 *
 * Multiplies the density and number of neighbours by the appropiate constants
 * and add the self-contribution term.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_end_volume(
    struct part* p) {

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float ih = 1.0f / h;

  /* Final operation on the density. */
  p->density.wcount =
      (p->density.wcount + kernel_root) * (4.0f / 3.0 * M_PI * kernel_gamma3);
  p->density.wcount_dh =
      p->density.wcount_dh * ih * (4.0f / 3.0 * M_PI * kernel_gamma3);
}

/**
 * @brief Prepare a particle for the force calculation.
 *
 * Computes viscosity term, conduction term and smoothing length gradient terms.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_gradient(
    struct part* p, struct xpart* xp) {

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float ih = 1.0f / h;
  const float ih2 = ih * ih;

  float detE, volume;
  float E[3][3];
  GFLOAT m, momentum[3], energy;

#ifndef THERMAL_ENERGY
  GFLOAT momentum2;
#endif

#if defined(SPH_GRADIENTS) && defined(SLOPE_LIMITER)
  GFLOAT gradrho[3], gradv[3][3], gradP[3];
  GFLOAT gradtrue, gradmax, gradmin, alpha;
#endif

  /* Final operation on the geometry. */
  /* we multiply with the smoothing kernel normalization ih3 and calculate the
   * volume */
  volume = ih * ih2 * (p->geometry.volume + kernel_root);
  p->geometry.volume = volume = 1. / volume;
  /* we multiply with the smoothing kernel normalization */
  p->geometry.matrix_E[0][0] = E[0][0] = ih * ih2 * p->geometry.matrix_E[0][0];
  p->geometry.matrix_E[0][1] = E[0][1] = ih * ih2 * p->geometry.matrix_E[0][1];
  p->geometry.matrix_E[0][2] = E[0][2] = ih * ih2 * p->geometry.matrix_E[0][2];
  p->geometry.matrix_E[1][0] = E[1][0] = ih * ih2 * p->geometry.matrix_E[1][0];
  p->geometry.matrix_E[1][1] = E[1][1] = ih * ih2 * p->geometry.matrix_E[1][1];
  p->geometry.matrix_E[1][2] = E[1][2] = ih * ih2 * p->geometry.matrix_E[1][2];
  p->geometry.matrix_E[2][0] = E[2][0] = ih * ih2 * p->geometry.matrix_E[2][0];
  p->geometry.matrix_E[2][1] = E[2][1] = ih * ih2 * p->geometry.matrix_E[2][1];
  p->geometry.matrix_E[2][2] = E[2][2] = ih * ih2 * p->geometry.matrix_E[2][2];

  /* invert the E-matrix */
  /* code shamelessly stolen from the public version of GIZMO */
  /* But since we should never invert a matrix, this code has to be replaced */
  detE = E[0][0] * E[1][1] * E[2][2] + E[0][1] * E[1][2] * E[2][0] +
         E[0][2] * E[1][0] * E[2][1] - E[0][2] * E[1][1] * E[2][0] -
         E[0][1] * E[1][0] * E[2][2] - E[0][0] * E[1][2] * E[2][1];
  /* check for zero determinant */
  if ((detE != 0) && !isnan(detE)) {
    p->geometry.matrix_E[0][0] = (E[1][1] * E[2][2] - E[1][2] * E[2][1]) / detE;
    p->geometry.matrix_E[0][1] = (E[0][2] * E[2][1] - E[0][1] * E[2][2]) / detE;
    p->geometry.matrix_E[0][2] = (E[0][1] * E[1][2] - E[0][2] * E[1][1]) / detE;
    p->geometry.matrix_E[1][0] = (E[1][2] * E[2][0] - E[1][0] * E[2][2]) / detE;
    p->geometry.matrix_E[1][1] = (E[0][0] * E[2][2] - E[0][2] * E[2][0]) / detE;
    p->geometry.matrix_E[1][2] = (E[0][2] * E[1][0] - E[0][0] * E[1][2]) / detE;
    p->geometry.matrix_E[2][0] = (E[1][0] * E[2][1] - E[1][1] * E[2][0]) / detE;
    p->geometry.matrix_E[2][1] = (E[0][1] * E[2][0] - E[0][0] * E[2][1]) / detE;
    p->geometry.matrix_E[2][2] = (E[0][0] * E[1][1] - E[0][1] * E[1][0]) / detE;
  } else {
    /* if the E-matrix is not well behaved, we cannot use it */
    p->geometry.matrix_E[0][0] = 0.0f;
    p->geometry.matrix_E[0][1] = 0.0f;
    p->geometry.matrix_E[0][2] = 0.0f;
    p->geometry.matrix_E[1][0] = 0.0f;
    p->geometry.matrix_E[1][1] = 0.0f;
    p->geometry.matrix_E[1][2] = 0.0f;
    p->geometry.matrix_E[2][0] = 0.0f;
    p->geometry.matrix_E[2][1] = 0.0f;
    p->geometry.matrix_E[2][2] = 0.0f;
  }

#ifdef SPH_GRADIENTS
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
#ifdef SLOPE_LIMITER
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
                   gradrho[2] * gradrho[2]);
  /* gradtrue might be zero. In this case, there is no gradient and we don't
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
#endif  // SLOPE_LIMITER
#else   // SPH_GRADIENTS
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
#endif  // SPH_GRADIENTS

  /* compute primitive variables */
  /* eqns (3)-(5) */
  m = p->conserved.mass;
  if (m) {
    momentum[0] = p->conserved.momentum[0];
    momentum[1] = p->conserved.momentum[1];
    momentum[2] = p->conserved.momentum[2];
#ifndef THERMAL_ENERGY
    momentum2 = (momentum[0] * momentum[0] + momentum[1] * momentum[1] +
                 momentum[2] * momentum[2]);
#endif
    energy = p->conserved.energy;
    p->primitives.rho = m / volume;
    p->primitives.v[0] = momentum[0] / m;
    p->primitives.v[1] = momentum[1] / m;
    p->primitives.v[2] = momentum[2] / m;
#ifndef THERMAL_ENERGY
    p->primitives.P =
        (const_hydro_gamma - 1.) * (energy - 0.5 * momentum2 / m) / volume;
#else
    p->primitives.P = (const_hydro_gamma - 1.) * energy / volume;
#endif
  }
}

/**
 * @brief Finishes the gradient calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_end_gradient(
    struct part* p) {

#ifndef SPH_GRADIENTS
  float h, ih, ih2, ih3;
#ifdef SLOPE_LIMITER
  GFLOAT gradrho[3], gradv[3][3], gradP[3];
  GFLOAT gradtrue, gradmax, gradmin, alpha;
#endif

  /* add kernel normalization to gradients */
  h = p->h;
  ih = 1.0f / h;
  ih2 = ih * ih;
  ih3 = ih * ih2;

  p->primitives.gradients.rho[0] *= ih3;
  p->primitives.gradients.rho[1] *= ih3;
  p->primitives.gradients.rho[2] *= ih3;

  p->primitives.gradients.v[0][0] *= ih3;
  p->primitives.gradients.v[0][1] *= ih3;
  p->primitives.gradients.v[0][2] *= ih3;
  p->primitives.gradients.v[1][0] *= ih3;
  p->primitives.gradients.v[1][1] *= ih3;
  p->primitives.gradients.v[1][2] *= ih3;
  p->primitives.gradients.v[2][0] *= ih3;
  p->primitives.gradients.v[2][1] *= ih3;
  p->primitives.gradients.v[2][2] *= ih3;

  p->primitives.gradients.P[0] *= ih3;
  p->primitives.gradients.P[1] *= ih3;
  p->primitives.gradients.P[2] *= ih3;

/* slope limiter */
#ifdef SLOPE_LIMITER
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

  gradtrue = gradrho[0] * gradrho[0] + gradrho[1] * gradrho[1] +
             gradrho[2] * gradrho[2];
  /* gradtrue might be zero. In this case, there is no gradient and we don't
     need to slope limit anything... */
  if (gradtrue) {
    gradtrue = sqrtf(gradtrue);
    gradtrue *= p->primitives.limiter.maxr;
    gradmax = p->primitives.limiter.rho[1] - p->primitives.rho;
    gradmin = p->primitives.rho - p->primitives.limiter.rho[0];
    /* gradmin and gradmax might be negative if the value of the current
       particle is larger/smaller than all neighbouring values */
    gradmax = fabs(gradmax);
    gradmin = fabs(gradmin);
    alpha = fmin(1.0f, fmin(gradmax / gradtrue, gradmin / gradtrue));
    p->primitives.gradients.rho[0] *= alpha;
    p->primitives.gradients.rho[1] *= alpha;
    p->primitives.gradients.rho[2] *= alpha;
  }

  gradtrue = gradv[0][0] * gradv[0][0] + gradv[0][1] * gradv[0][1] +
             gradv[0][2] * gradv[0][2];
  if (gradtrue) {
    gradtrue = sqrtf(gradtrue);
    gradtrue *= p->primitives.limiter.maxr;
    gradmax = p->primitives.limiter.v[0][1] - p->primitives.v[0];
    gradmin = p->primitives.v[0] - p->primitives.limiter.v[0][0];
    gradmax = fabs(gradmax);
    gradmin = fabs(gradmin);
    alpha = fmin(1.0f, fmin(gradmax / gradtrue, gradmin / gradtrue));
    p->primitives.gradients.v[0][0] *= alpha;
    p->primitives.gradients.v[0][1] *= alpha;
    p->primitives.gradients.v[0][2] *= alpha;
  }

  gradtrue = gradv[1][0] * gradv[1][0] + gradv[1][1] * gradv[1][1] +
             gradv[1][2] * gradv[1][2];
  if (gradtrue) {
    gradtrue = sqrtf(gradtrue);
    gradtrue *= p->primitives.limiter.maxr;
    gradmax = p->primitives.limiter.v[1][1] - p->primitives.v[1];
    gradmin = p->primitives.v[1] - p->primitives.limiter.v[1][0];
    gradmax = fabs(gradmax);
    gradmin = fabs(gradmin);
    alpha = fmin(1.0f, fmin(gradmax / gradtrue, gradmin / gradtrue));
    p->primitives.gradients.v[1][0] *= alpha;
    p->primitives.gradients.v[1][1] *= alpha;
    p->primitives.gradients.v[1][2] *= alpha;
  }

  gradtrue = gradv[2][0] * gradv[2][0] + gradv[2][1] * gradv[2][1] +
             gradv[2][2] * gradv[2][2];
  if (gradtrue) {
    gradtrue = sqrtf(gradtrue);
    gradtrue *= p->primitives.limiter.maxr;
    gradmax = p->primitives.limiter.v[2][1] - p->primitives.v[2];
    gradmin = p->primitives.v[2] - p->primitives.limiter.v[2][0];
    gradmax = fabs(gradmax);
    gradmin = fabs(gradmin);
    alpha = fmin(1.0f, fmin(gradmax / gradtrue, gradmin / gradtrue));
    p->primitives.gradients.v[2][0] *= alpha;
    p->primitives.gradients.v[2][1] *= alpha;
    p->primitives.gradients.v[2][2] *= alpha;
  }

  gradtrue = gradP[0] * gradP[0] + gradP[1] * gradP[1] + gradP[2] * gradP[2];
  if (gradtrue) {
    gradtrue = sqrtf(gradtrue);
    gradtrue *= p->primitives.limiter.maxr;
    gradmax = p->primitives.limiter.P[1] - p->primitives.P;
    gradmin = p->primitives.P - p->primitives.limiter.P[0];
    gradmax = fabs(gradmax);
    gradmin = fabs(gradmin);
    alpha = fmin(1.0f, fmin(gradmax / gradtrue, gradmin / gradtrue));
    p->primitives.gradients.P[0] *= alpha;
    p->primitives.gradients.P[1] *= alpha;
    p->primitives.gradients.P[2] *= alpha;
  }
#endif  // SLOPE_LIMITER

#endif  // SPH_GRADIENTS
}

/**
 * @brief Prepare a particle for the fluxes calculation.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_fluxes(
    struct part* p, struct xpart* xp) {

  /* initialize variables used for timestep calculation */
  p->timestepvars.vmax = 0.0f;
}

/**
 * @brief Reset acceleration fields of a particle
 *
 * Resets all hydro acceleration and time derivative fields in preparation
 * for the sums taking place in the variaous force tasks
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_reset_acceleration(
    struct part* p) {

  /* figure out what to put here */
}

/**
 * @brief Finishes the fluxes calculation.
 *
 * Multiplies the forces and accelerationsby the appropiate constants
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_end_fluxes(
    struct part* p) {

  /* do nothing */
}

/**
 * @brief Converts hydro quantity of a particle
 *
 * Requires the volume to be known
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_convert_quantities(
    struct part* p) {

  float volume;
  GFLOAT m;
  GFLOAT momentum[3];
#ifndef THERMAL_ENERGY
  GFLOAT momentum2;
#endif
  volume = p->geometry.volume;

  /* set hydro velocities */
  p->primitives.v[0] = p->v[0];
  p->primitives.v[1] = p->v[1];
  p->primitives.v[2] = p->v[2];
  /* P actually contains internal energy at this point */
  p->primitives.P *= (const_hydro_gamma - 1.) * p->primitives.rho;

  p->conserved.mass = m = p->primitives.rho * volume;
  p->conserved.momentum[0] = momentum[0] = m * p->primitives.v[0];
  p->conserved.momentum[1] = momentum[1] = m * p->primitives.v[1];
  p->conserved.momentum[2] = momentum[2] = m * p->primitives.v[2];
#ifndef THERMAL_ENERGY
  momentum2 = momentum[0] * momentum[0] + momentum[1] * momentum[1] +
              momentum[2] * momentum[2];
  p->conserved.energy =
      p->primitives.P / (const_hydro_gamma - 1.) * volume + 0.5 * momentum2 / m;
#else
  p->conserved.energy = p->primitives.P / (const_hydro_gamma - 1.) * volume;
#endif
}

// MATTHIEU
__attribute__((always_inline)) INLINE static void hydro_end_density(
    struct part* p, float time) {}
__attribute__((always_inline)) INLINE static void hydro_prepare_force(
    struct part* p, struct xpart* xp, int ti_current, double timeBase) {}
__attribute__((always_inline)) INLINE static void hydro_predict_extra(
    struct part* p, struct xpart* xp, int t0, int t1, double timeBase) {}
__attribute__((always_inline)) INLINE static void hydro_end_force(
    struct part* p) {}
__attribute__((always_inline)) INLINE static void hydro_kick_extra(
    struct part* p, struct xpart* xp, float dt, float half_dt) {}
__attribute__((always_inline)) INLINE static float hydro_get_internal_energy(
    struct part* p) {
  return 0.f;
}
