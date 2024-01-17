/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_GIZMO_HYDRO_IACT_H
#define SWIFT_GIZMO_HYDRO_IACT_H

#include "chemistry_additions.h"
#include "hydro_flux.h"
#include "hydro_getters.h"
#include "hydro_gradients.h"
#include "hydro_setters.h"
#include "hydro_velocities.h"
#include "rt_additions.h"

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
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  float wi, wj, wi_dx, wj_dx;

  /* Get r and h inverse. */
  const float r = sqrtf(r2);

  /* Compute density of pi. */
  const float hi_inv = 1.0f / hi;
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + xi * wi_dx);

  /* these are eqns. (1) and (2) in the summary */
  pi->geometry.volume += wi;
  for (int k = 0; k < 3; k++)
    for (int l = 0; l < 3; l++)
      pi->geometry.matrix_E[k][l] += dx[k] * dx[l] * wi;

  hydro_velocities_update_centroid_left(pi, dx, wi);

  /* Compute density of pj. */
  const float hj_inv = 1.0f / hj;
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);

  pj->density.wcount += wj;
  pj->density.wcount_dh -= (hydro_dimension * wj + xj * wj_dx);

  /* these are eqns. (1) and (2) in the summary */
  pj->geometry.volume += wj;
  for (int k = 0; k < 3; k++)
    for (int l = 0; l < 3; l++)
      pj->geometry.matrix_E[k][l] += dx[k] * dx[l] * wj;

  hydro_velocities_update_centroid_right(pj, dx, wj);
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
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj, const float a,
    const float H) {

  float wi, wi_dx;

  /* Get r and h inverse. */
  const float r = sqrtf(r2);

  const float hi_inv = 1.0f / hi;
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + xi * wi_dx);

  /* these are eqns. (1) and (2) in the summary */
  pi->geometry.volume += wi;
  for (int k = 0; k < 3; k++)
    for (int l = 0; l < 3; l++)
      pi->geometry.matrix_E[k][l] += dx[k] * dx[l] * wi;

  hydro_velocities_update_centroid_left(pi, dx, wi);
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j
 *
 * This method wraps around hydro_gradients_collect, which can be an empty
 * method, in which case no gradients are used.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  hydro_gradients_collect(r2, dx, hi, hj, pi, pj);
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j:
 * non-symmetric version
 *
 * This method wraps around hydro_gradients_nonsym_collect, which can be an
 * empty method, in which case no gradients are used.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  hydro_gradients_nonsym_collect(r2, dx, hi, hj, pi, pj);
}

/**
 * @brief Common part of the flux calculation between particle i and j
 *
 * Since the only difference between the symmetric and non-symmetric version
 * of the flux calculation  is in the update of the conserved variables at the
 * very end (which is not done for particle j if mode is 0), both
 * runner_iact_force and runner_iact_nonsym_force call this method, with an
 * appropriate mode.
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
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_fluxes_common(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, int mode, const float a,
    const float H) {

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Initialize local variables */
  float Bi[3][3];
  float Bj[3][3];
  float vi[3], vj[3];
  for (int k = 0; k < 3; k++) {
    for (int l = 0; l < 3; l++) {
      Bi[k][l] = pi->geometry.matrix_E[k][l];
      Bj[k][l] = pj->geometry.matrix_E[k][l];
    }
    vi[k] = pi->v[k]; /* particle velocities */
    vj[k] = pj->v[k];
  }
  const float Vi = pi->geometry.volume;
  const float Vj = pj->geometry.volume;
  float Wi[5], Wj[5];
  hydro_part_get_primitive_variables(pi, Wi);
  hydro_part_get_primitive_variables(pj, Wj);

  /* calculate the maximal signal velocity */
  float vmax;
  if (Wi[0] > 0.0f && Wj[0] > 0.0f) {
    const float ci = gas_soundspeed_from_pressure(Wi[0], Wi[4]);
    const float cj = gas_soundspeed_from_pressure(Wj[0], Wj[4]);
    vmax = ci + cj;
  } else {
    vmax = 0.0f;
  }

  float dvdr = (pi->v[0] - pj->v[0]) * dx[0] + (pi->v[1] - pj->v[1]) * dx[1] +
               (pi->v[2] - pj->v[2]) * dx[2];

  /* Velocity on the axis linking the particles */
  /* This velocity will be the same as dvdr for MFM, so hopefully this gets
     optimised out. */
  float dvdotdx = (Wi[1] - Wj[1]) * dx[0] + (Wi[2] - Wj[2]) * dx[1] +
                  (Wi[3] - Wj[3]) * dx[2];

  /* We only care about this velocity for particles moving towards each others
   */
  dvdotdx = min3(dvdr, dvdotdx, 0.f);

  /* Get the signal velocity */
  vmax -= const_viscosity_beta * dvdotdx * r_inv;

  /* Store the signal velocity */
  pi->timestepvars.vmax = max(pi->timestepvars.vmax, vmax);
  if (mode == 1) {
    pj->timestepvars.vmax = max(pj->timestepvars.vmax, vmax);
  }

  /* Compute kernel of pi. */
  float wi, wi_dx;
  const float hi_inv = 1.0f / hi;
  const float hi_inv_dim = pow_dimension(hi_inv);
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  /* Compute kernel of pj. */
  float wj, wj_dx;
  const float hj_inv = 1.0f / hj;
  const float hj_inv_dim = pow_dimension(hj_inv);
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);

  /* Compute h_dt. We are going to use an SPH-like estimate of div_v for that */
  const float hidp1 = pow_dimension_plus_one(hi_inv);
  const float hjdp1 = pow_dimension_plus_one(hj_inv);
  const float wi_dr = hidp1 * wi_dx;
  const float wj_dr = hjdp1 * wj_dx;
  dvdr *= r_inv;
  if (Wj[0] > 0.0f) {
    pi->force.h_dt -= pj->conserved.mass * dvdr / Wj[0] * wi_dr;
  }
  if (mode == 1 && Wi[0] > 0.0f) {
    pj->force.h_dt -= pi->conserved.mass * dvdr / Wi[0] * wj_dr;
  }

  /* Compute (square of) area */
  /* eqn. (7) */
  float Anorm2 = 0.0f;
  float A[3];
  if (hydro_part_geometry_well_behaved(pi) &&
      hydro_part_geometry_well_behaved(pj)) {
    /* in principle, we use Vi and Vj as weights for the left and right
       contributions to the generalized surface vector.
       However, if Vi and Vj are very different (because they have very
       different
       smoothing lengths), then the expressions below are more stable. */
    float Xi = Vi;
    float Xj = Vj;
#ifdef GIZMO_VOLUME_CORRECTION
    if (fabsf(Vi - Vj) / min(Vi, Vj) > 1.5f * hydro_dimension) {
      Xi = (Vi * hj + Vj * hi) / (hi + hj);
      Xj = Xi;
    }
#endif
    for (int k = 0; k < 3; k++) {
      /* we add a minus sign since dx is pi->x - pj->x */
      A[k] = -Xi * (Bi[k][0] * dx[0] + Bi[k][1] * dx[1] + Bi[k][2] * dx[2]) *
                 wi * hi_inv_dim -
             Xj * (Bj[k][0] * dx[0] + Bj[k][1] * dx[1] + Bj[k][2] * dx[2]) *
                 wj * hj_inv_dim;
      Anorm2 += A[k] * A[k];
    }
  } else {
    /* ill condition gradient matrix: revert to SPH face area */
    const float Anorm =
        -(hidp1 * Vi * Vi * wi_dx + hjdp1 * Vj * Vj * wj_dx) * r_inv;
    A[0] = -Anorm * dx[0];
    A[1] = -Anorm * dx[1];
    A[2] = -Anorm * dx[2];
    Anorm2 = Anorm * Anorm * r2;
  }

  /* if the interface has no area, nothing happens and we return */
  /* continuing results in dividing by zero and NaN's... */
  if (Anorm2 == 0.0f) {
    return;
  }

  /* Compute the area */
  const float Anorm_inv = 1.0f / sqrtf(Anorm2);
  const float Anorm = Anorm2 * Anorm_inv;

#ifdef SWIFT_DEBUG_CHECKS
  /* For stability reasons, we do require A and dx to have opposite
     directions (basically meaning that the surface normal for the surface
     always points from particle i to particle j, as it would in a real
     moving-mesh code). If not, our scheme is no longer upwind and hence can
     become unstable. */
  const float dA_dot_dx = A[0] * dx[0] + A[1] * dx[1] + A[2] * dx[2];
  /* In GIZMO, Phil Hopkins reverts to an SPH integration scheme if this
     happens. We curently just ignore this case and display a message. */
  const float rdim = pow_dimension(r);
  if (dA_dot_dx > 1.e-6f * rdim) {
    message("Ill conditioned gradient matrix (%g %g %g %g %g)!", dA_dot_dx,
            Anorm, Vi, Vj, r);
  }
#endif

  /* compute the normal vector of the interface */
  const float n_unit[3] = {A[0] * Anorm_inv, A[1] * Anorm_inv,
                           A[2] * Anorm_inv};

  /* Compute interface position (relative to pi, since we don't need the actual
   * position) eqn. (8) */
  const float xfac = -hi / (hi + hj);
  const float xij_i[3] = {xfac * dx[0], xfac * dx[1], xfac * dx[2]};

  /* Compute interface velocity */
  /* eqn. (9) */
  const float vij[3] = {vi[0] + (vi[0] - vj[0]) * xfac,
                        vi[1] + (vi[1] - vj[1]) * xfac,
                        vi[2] + (vi[2] - vj[2]) * xfac};

  /* complete calculation of position of interface */
  /* NOTE: dx is not necessarily just pi->x - pj->x but can also contain
           correction terms for periodicity. If we do the interpolation,
           we have to use xij w.r.t. the actual particle.
           => we need a separate xij for pi and pj... */
  /* tldr: we do not need the code below, but we do need the same code as above
     but then with i and j swapped */
  //    for ( k = 0 ; k < 3 ; k++ )
  //      xij[k] += pi->x[k];

  hydro_gradients_predict(pi, pj, hi, hj, dx, r, xij_i, Wi, Wj);

  /* Boost the primitive variables to the frame of reference of the interface */
  /* Note that velocities are indices 1-3 in W */
  Wi[1] -= vij[0];
  Wi[2] -= vij[1];
  Wi[3] -= vij[2];
  Wj[1] -= vij[0];
  Wj[2] -= vij[1];
  Wj[3] -= vij[2];

  /* we don't need to rotate, we can use the unit vector in the Riemann problem
   * itself (see GIZMO) */

  float totflux[5];
  hydro_compute_flux(Wi, Wj, n_unit, vij, Anorm, totflux);

  /* get the time step for the flux exchange. This is always the smallest time
     step among the two particles */
  const float mindt =
      (pj->flux.dt > 0.0f) ? fminf(pi->flux.dt, pj->flux.dt) : pi->flux.dt;

  hydro_part_update_fluxes_left(pi, totflux, dx, mindt);

  /* Unlike in SPH schemes, we do need to update inactive neighbours, so that
     the fluxes are always exchanged symmetrically. Thanks to our sneaky use
     of flux.dt, we can detect inactive neighbours through their negative
     time step. */
  if (mode == 1 || (pj->flux.dt < 0.0f)) {
    hydro_part_update_fluxes_right(pj, totflux, dx, mindt);
  }

  /* If we're working with RT, we need to pay additional attention to the
   * individual mass fractions of ionizing species. */
  rt_part_update_mass_fluxes(pi, pj, totflux[0], mode);

  /* Advect metals if working with chemistry. */
  runner_iact_chemistry_fluxes(pi, pj, totflux[0], mindt, mode);
}

/**
 * @brief Flux calculation between particle i and particle j
 *
 * This method calls runner_iact_fluxes_common with mode 1.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  runner_iact_fluxes_common(r2, dx, hi, hj, pi, pj, 1, a, H);
}

/**
 * @brief Flux calculation between particle i and particle j: non-symmetric
 * version
 *
 * This method calls runner_iact_fluxes_common with mode 0.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  runner_iact_fluxes_common(r2, dx, hi, hj, pi, pj, 0, a, H);
}

#endif /* SWIFT_GIZMO_HYDRO_IACT_H */
