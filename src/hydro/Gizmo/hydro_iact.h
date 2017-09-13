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

#include "adiabatic_index.h"
#include "hydro_flux_limiters.h"
#include "hydro_gradients.h"
#include "riemann.h"

#define GIZMO_VOLUME_CORRECTION

/**
 * @brief Calculate the volume interaction between particle i and particle j
 *
 * The volume is in essence the same as the weighted number of neighbours in a
 * classical SPH density calculation.
 *
 * We also calculate the components of the matrix E, which is used for second
 * order accurate gradient calculations and for the calculation of the interface
 * surface areas.
 *
 * @param r2 Squared distance between particle i and particle j.
 * @param dx Distance vector between the particles (dx = pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
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
  pi->density.wcount_dh -= (hydro_dimension * wi + xi * wi_dx);

  /* these are eqns. (1) and (2) in the summary */
  pi->geometry.volume += wi;
  for (k = 0; k < 3; k++)
    for (l = 0; l < 3; l++) pi->geometry.matrix_E[k][l] += dx[k] * dx[l] * wi;

  pi->geometry.centroid[0] -= dx[0] * wi;
  pi->geometry.centroid[1] -= dx[1] * wi;
  pi->geometry.centroid[2] -= dx[2] * wi;

  /* Compute density of pj. */
  h_inv = 1.0 / hj;
  xj = r * h_inv;
  kernel_deval(xj, &wj, &wj_dx);

  pj->density.wcount += wj;
  pj->density.wcount_dh -= (hydro_dimension * wj + xj * wj_dx);

  /* these are eqns. (1) and (2) in the summary */
  pj->geometry.volume += wj;
  for (k = 0; k < 3; k++)
    for (l = 0; l < 3; l++) pj->geometry.matrix_E[k][l] += dx[k] * dx[l] * wj;

  pj->geometry.centroid[0] += dx[0] * wj;
  pj->geometry.centroid[1] += dx[1] * wj;
  pj->geometry.centroid[2] += dx[2] * wj;
}

/**
 * @brief Calculate the volume interaction between particle i and particle j:
 * non-symmetric version
 *
 * The volume is in essence the same as the weighted number of neighbours in a
 * classical SPH density calculation.
 *
 * We also calculate the components of the matrix E, which is used for second
 * order accurate gradient calculations and for the calculation of the interface
 * surface areas.
 *
 * @param r2 Squared distance between particle i and particle j.
 * @param dx Distance vector between the particles (dx = pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

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
  pi->density.wcount_dh -= (hydro_dimension * wi + xi * wi_dx);

  /* these are eqns. (1) and (2) in the summary */
  pi->geometry.volume += wi;
  for (k = 0; k < 3; k++)
    for (l = 0; l < 3; l++) pi->geometry.matrix_E[k][l] += dx[k] * dx[l] * wi;

  pi->geometry.centroid[0] -= dx[0] * wi;
  pi->geometry.centroid[1] -= dx[1] * wi;
  pi->geometry.centroid[2] -= dx[2] * wi;
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j
 *
 * This method wraps around hydro_gradients_collect, which can be an empty
 * method, in which case no gradients are used.
 *
 * @param r2 Squared distance between particle i and particle j.
 * @param dx Distance vector between the particles (dx = pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void runner_iact_gradient(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float hi_inv, hi_inv_dim, xi, wi, wi_dx;
  float hj_inv, hj_inv_dim, xj, wj, wj_dx;
  float Bi[3][3], Bj[3][3];
  float Vi, Vj;
  float A, Anorm;
  int k, l;
  float r;

  r = sqrtf(r2);

  hi_inv = 1.0 / hi;
  hi_inv_dim = pow_dimension(hi_inv);
  xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  /* Compute kernel of pj. */
  hj_inv = 1.0 / hj;
  hj_inv_dim = pow_dimension(hj_inv);
  xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);

  for (k = 0; k < 3; k++) {
    for (l = 0; l < 3; l++) {
      Bi[k][l] = pi->geometry.matrix_E[k][l];
      Bj[k][l] = pj->geometry.matrix_E[k][l];
    }
  }
  Vi = pi->geometry.volume;
  Vj = pj->geometry.volume;

  /* Compute area */
  /* eqn. (7) */
  Anorm = 0.0f;
  for (k = 0; k < 3; k++) {
    /* we add a minus sign since dx is pi->x - pj->x */
    A = -Vi * (Bi[k][0] * dx[0] + Bi[k][1] * dx[1] + Bi[k][2] * dx[2]) * wi *
            hi_inv_dim -
        Vj * (Bj[k][0] * dx[0] + Bj[k][1] * dx[1] + Bj[k][2] * dx[2]) * wj *
            hj_inv_dim;
    Anorm += A * A;
  }

  Anorm = sqrtf(Anorm);

  pi->geometry.Atot += Anorm;
  pj->geometry.Atot += Anorm;

  hydro_gradients_collect(r2, dx, hi, hj, pi, pj);
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j:
 * non-symmetric version
 *
 * This method wraps around hydro_gradients_nonsym_collect, which can be an
 * empty method, in which case no gradients are used.
 *
 * @param r2 Squared distance between particle i and particle j.
 * @param dx Distance vector between the particles (dx = pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_gradient(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float hi_inv, hi_inv_dim, xi, wi, wi_dx;
  float hj_inv, hj_inv_dim, xj, wj, wj_dx;
  float Bi[3][3], Bj[3][3];
  float Vi, Vj;
  float A, Anorm;
  int k, l;
  float r;

  r = sqrtf(r2);

  hi_inv = 1.0 / hi;
  hi_inv_dim = pow_dimension(hi_inv);
  xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  /* Compute kernel of pj. */
  hj_inv = 1.0 / hj;
  hj_inv_dim = pow_dimension(hj_inv);
  xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);

  for (k = 0; k < 3; k++) {
    for (l = 0; l < 3; l++) {
      Bi[k][l] = pi->geometry.matrix_E[k][l];
      Bj[k][l] = pj->geometry.matrix_E[k][l];
    }
  }
  Vi = pi->geometry.volume;
  Vj = pj->geometry.volume;

  /* Compute area */
  /* eqn. (7) */
  Anorm = 0.0f;
  for (k = 0; k < 3; k++) {
    /* we add a minus sign since dx is pi->x - pj->x */
    A = -Vi * (Bi[k][0] * dx[0] + Bi[k][1] * dx[1] + Bi[k][2] * dx[2]) * wi *
            hi_inv_dim -
        Vj * (Bj[k][0] * dx[0] + Bj[k][1] * dx[1] + Bj[k][2] * dx[2]) * wj *
            hj_inv_dim;
    Anorm += A * A;
  }

  Anorm = sqrtf(Anorm);

  pi->geometry.Atot += Anorm;

  hydro_gradients_nonsym_collect(r2, dx, hi, hj, pi, pj);
}

/**
 * @brief Common part of the flux calculation between particle i and j
 *
 * Since the only difference between the symmetric and non-symmetric version
 * of the flux calculation  is in the update of the conserved variables at the
 * very end (which is not done for particle j if mode is 0 and particle j is
 * active), both runner_iact_force and runner_iact_nonsym_force call this
 * method, with an appropriate mode.
 *
 * This method calculates the surface area of the interface between particle i
 * and particle j, as well as the interface position and velocity. These are
 * then used to reconstruct and predict the primitive variables, which are then
 * fed to a Riemann solver that calculates a flux. This flux is used to update
 * the conserved variables of particle i or both particles.
 *
 * This method also calculates the maximal velocity used to calculate the time
 * step.
 *
 * @param r2 Squared distance between particle i and particle j.
 * @param dx Distance vector between the particles (dx = pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void runner_iact_fluxes_common(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj,
    int mode) {

  float r = sqrtf(r2);
  float xi, xj;
  float hi_inv, hi_inv_dim;
  float hj_inv, hj_inv_dim;
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
  float Wi[5], Wj[5];
  float dti, dtj, mindt;
  float n_unit[3];

  /* Initialize local variables */
  for (k = 0; k < 3; k++) {
    for (l = 0; l < 3; l++) {
      Bi[k][l] = pi->geometry.matrix_E[k][l];
      Bj[k][l] = pj->geometry.matrix_E[k][l];
    }
    vi[k] = pi->force.v_full[k]; /* particle velocities */
    vj[k] = pj->force.v_full[k];
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

  dti = pi->force.dt;
  dtj = pj->force.dt;

  /* calculate the maximal signal velocity */
  if (Wi[0] > 0.0f && Wj[0] > 0.0f) {
    vmax =
        sqrtf(hydro_gamma * Wi[4] / Wi[0]) + sqrtf(hydro_gamma * Wj[4] / Wj[0]);
  } else {
    vmax = 0.0f;
  }
  dvdotdx = (Wi[1] - Wj[1]) * dx[0] + (Wi[2] - Wj[2]) * dx[1] +
            (Wi[3] - Wj[3]) * dx[2];
  dvdotdx = min(dvdotdx, (vi[0] - vj[0]) * dx[0] + (vi[1] - vj[1]) * dx[1] +
                             (vi[2] - vj[2]) * dx[2]);
  if (dvdotdx < 0.) {
    /* the magical factor 3 also appears in Gadget2 */
    vmax -= 3. * dvdotdx / r;
  }
  pi->timestepvars.vmax = max(pi->timestepvars.vmax, vmax);
  if (mode == 1) {
    pj->timestepvars.vmax = max(pj->timestepvars.vmax, vmax);
  }

  /* The flux will be exchanged using the smallest time step of the two
   * particles */
  mindt = min(dti, dtj);

  /* Compute kernel of pi. */
  hi_inv = 1.0 / hi;
  hi_inv_dim = pow_dimension(hi_inv);
  xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  /* Compute kernel of pj. */
  hj_inv = 1.0 / hj;
  hj_inv_dim = pow_dimension(hj_inv);
  xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);

  /* Compute h_dt. We are going to use an SPH-like estimate of div_v for that */
  float dvdr = (pi->v[0] - pj->v[0]) * dx[0] + (pi->v[1] - pj->v[1]) * dx[1] +
               (pi->v[2] - pj->v[2]) * dx[2];
  float ri = 1.0f / r;
  float hidp1 = pow_dimension_plus_one(hi_inv);
  float hjdp1 = pow_dimension_plus_one(hj_inv);
  float wi_dr = hidp1 * wi_dx;
  float wj_dr = hjdp1 * wj_dx;
  dvdr *= ri;
  if (pj->primitives.rho > 0.) {
    pi->force.h_dt -= pj->conserved.mass * dvdr / pj->primitives.rho * wi_dr;
  }
  if (mode == 1 && pi->primitives.rho > 0.) {
    pj->force.h_dt -= pi->conserved.mass * dvdr / pi->primitives.rho * wj_dr;
  }

  /* Compute area */
  /* eqn. (7) */
  Anorm = 0.0f;
  if (pi->density.wcorr > const_gizmo_min_wcorr &&
      pj->density.wcorr > const_gizmo_min_wcorr) {
    /* in principle, we use Vi and Vj as weights for the left and right
       contributions to the generalized surface vector.
       However, if Vi and Vj are very different (because they have very
       different
       smoothing lengths), then the expressions below are more stable. */
    float Xi = Vi;
    float Xj = Vj;
#ifdef GIZMO_VOLUME_CORRECTION
    if (fabsf(Vi - Vj) / fminf(Vi, Vj) > 1.5 * hydro_dimension) {
      Xi = (Vi * hj + Vj * hi) / (hi + hj);
      Xj = Xi;
    }
#endif
    for (k = 0; k < 3; k++) {
      /* we add a minus sign since dx is pi->x - pj->x */
      A[k] = -Xi * (Bi[k][0] * dx[0] + Bi[k][1] * dx[1] + Bi[k][2] * dx[2]) *
                 wj * hj_inv_dim -
             Xj * (Bj[k][0] * dx[0] + Bj[k][1] * dx[1] + Bj[k][2] * dx[2]) *
                 wi * hi_inv_dim;
      Anorm += A[k] * A[k];
    }
  } else {
    /* ill condition gradient matrix: revert to SPH face area */
    Anorm = -(hidp1 * Vi * Vi * wi_dx + hjdp1 * Vj * Vj * wj_dx) * ri;
    A[0] = -Anorm * dx[0];
    A[1] = -Anorm * dx[1];
    A[2] = -Anorm * dx[2];
    Anorm *= Anorm * r2;
  }

  if (Anorm == 0.) {
    /* if the interface has no area, nothing happens and we return */
    /* continuing results in dividing by zero and NaN's... */
    return;
  }

  Anorm = sqrtf(Anorm);

#ifdef SWIFT_DEBUG_CHECKS
  /* For stability reasons, we do require A and dx to have opposite
     directions (basically meaning that the surface normal for the surface
     always points from particle i to particle j, as it would in a real
     moving-mesh code). If not, our scheme is no longer upwind and hence can
     become unstable. */
  float dA_dot_dx = A[0] * dx[0] + A[1] * dx[1] + A[2] * dx[2];
  /* In GIZMO, Phil Hopkins reverts to an SPH integration scheme if this
     happens. We curently just ignore this case and display a message. */
  const float rdim = pow_dimension(r);
  if (dA_dot_dx > 1.e-6 * rdim) {
    message("Ill conditioned gradient matrix (%g %g %g %g %g)!", dA_dot_dx,
            Anorm, Vi, Vj, r);
  }
#endif

  /* compute the normal vector of the interface */
  for (k = 0; k < 3; k++) n_unit[k] = A[k] / Anorm;

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

  hydro_gradients_predict(pi, pj, hi, hj, dx, r, xij_i, Wi, Wj, mindt);

  /* we don't need to rotate, we can use the unit vector in the Riemann problem
   * itself (see GIZMO) */

  float totflux[5];
  riemann_solve_for_flux(Wi, Wj, n_unit, vij, totflux);

  hydro_flux_limiters_apply(totflux, pi, pj);

  /* Store mass flux */
  float mflux = Anorm * totflux[0];
  pi->gravity.mflux[0] += mflux * dx[0];
  pi->gravity.mflux[1] += mflux * dx[1];
  pi->gravity.mflux[2] += mflux * dx[2];

  /* Update conserved variables */
  /* eqn. (16) */
  pi->conserved.flux.mass -= mindt * Anorm * totflux[0];
  pi->conserved.flux.momentum[0] -= mindt * Anorm * totflux[1];
  pi->conserved.flux.momentum[1] -= mindt * Anorm * totflux[2];
  pi->conserved.flux.momentum[2] -= mindt * Anorm * totflux[3];
  pi->conserved.flux.energy -= mindt * Anorm * totflux[4];

#ifndef GIZMO_TOTAL_ENERGY
  float ekin = 0.5f * (pi->primitives.v[0] * pi->primitives.v[0] +
                       pi->primitives.v[1] * pi->primitives.v[1] +
                       pi->primitives.v[2] * pi->primitives.v[2]);
  pi->conserved.flux.energy += mindt * Anorm * totflux[1] * pi->primitives.v[0];
  pi->conserved.flux.energy += mindt * Anorm * totflux[2] * pi->primitives.v[1];
  pi->conserved.flux.energy += mindt * Anorm * totflux[3] * pi->primitives.v[2];
  pi->conserved.flux.energy -= mindt * Anorm * totflux[0] * ekin;
#endif

  /* here is how it works:
     Mode will only be 1 if both particles are ACTIVE and they are in the same
     cell. In this case, this method IS the flux calculation for particle j, and
     we HAVE TO UPDATE it.
     Mode 0 can mean several things: it can mean that particle j is INACTIVE, in
     which case we NEED TO UPDATE it, since otherwise the flux is lost from the
     system and the conserved variable is not conserved.
     It can also mean that particle j sits in another cell and is ACTIVE. In
     this case, the flux exchange for particle j is done TWICE and we SHOULD NOT
     UPDATE particle j.
     ==> we update particle j if (MODE IS 1) OR (j IS INACTIVE)
  */

  if (mode == 1 || pj->force.active == 0) {
    /* Store mass flux */
    mflux = Anorm * totflux[0];
    pj->gravity.mflux[0] -= mflux * dx[0];
    pj->gravity.mflux[1] -= mflux * dx[1];
    pj->gravity.mflux[2] -= mflux * dx[2];

    pj->conserved.flux.mass += mindt * Anorm * totflux[0];
    pj->conserved.flux.momentum[0] += mindt * Anorm * totflux[1];
    pj->conserved.flux.momentum[1] += mindt * Anorm * totflux[2];
    pj->conserved.flux.momentum[2] += mindt * Anorm * totflux[3];
    pj->conserved.flux.energy += mindt * Anorm * totflux[4];

#ifndef GIZMO_TOTAL_ENERGY
    ekin = 0.5f * (pj->primitives.v[0] * pj->primitives.v[0] +
                   pj->primitives.v[1] * pj->primitives.v[1] +
                   pj->primitives.v[2] * pj->primitives.v[2]);
    pj->conserved.flux.energy -=
        mindt * Anorm * totflux[1] * pj->primitives.v[0];
    pj->conserved.flux.energy -=
        mindt * Anorm * totflux[2] * pj->primitives.v[1];
    pj->conserved.flux.energy -=
        mindt * Anorm * totflux[3] * pj->primitives.v[2];
    pj->conserved.flux.energy += mindt * Anorm * totflux[0] * ekin;
#endif
  }
}

/**
 * @brief Flux calculation between particle i and particle j
 *
 * This method calls runner_iact_fluxes_common with mode 1.
 *
 * @param r2 Squared distance between particle i and particle j.
 * @param dx Distance vector between the particles (dx = pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  runner_iact_fluxes_common(r2, dx, hi, hj, pi, pj, 1);
}

/**
 * @brief Flux calculation between particle i and particle j: non-symmetric
 * version
 *
 * This method calls runner_iact_fluxes_common with mode 0.
 *
 * @param r2 Squared distance between particle i and particle j.
 * @param dx Distance vector between the particles (dx = pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  runner_iact_fluxes_common(r2, dx, hi, hj, pi, pj, 0);
}
